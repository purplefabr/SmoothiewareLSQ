#include "LSQDeltaCalibrationStrategy.h"
#include "Kernel.h"
#include "Config.h"
#include "Robot.h"
#include "StreamOutputPool.h"
#include "Gcode.h"
#include "checksumm.h"
#include "ConfigValue.h"
#include "PublicDataRequest.h"
#include "EndstopsPublicAccess.h"
#include "PublicData.h"
#include "Conveyor.h"
#include "ZProbe.h"
#include "BaseSolution.h"
#include "Matrix.h"
#include "ActuatorCoordinates.h"

#include <fastmath.h>
#include <tuple>
#include <vector>

#define PIOVER180   0.01745329251994329576923690768489F
#define PI 3.14159265358979F

#define sample_count_checksum CHECKSUM("sample_count")
#define factors_checksum CHECKSUM("factors")
#define radius_checksum       CHECKSUM("radius")
#define initial_height_checksum CHECKSUM("initial_height")

// deprecated
#define probe_radius_checksum CHECKSUM("probe_radius")

const int MAX_DELTA_PROBE_POINTS = 13;
const int NumDeltaFactors = 7;

bool LSQDeltaCalibrationStrategy::handleConfig()
{
    // default is probably wrong
    float r = THEKERNEL->config->value(leveling_strategy_checksum, lsq_delta_calibration_strategy_checksum, radius_checksum)->by_default(-1)->as_number();
    if (r == -1) {
        // deprecated config syntax]
        r = THEKERNEL->config->value(zprobe_checksum, probe_radius_checksum)->by_default(100.0F)->as_number();
    }
    this->probe_radius = r;

    this->sample_count = THEKERNEL->config->value(leveling_strategy_checksum, lsq_delta_calibration_strategy_checksum, sample_count_checksum)->by_default(13)->as_number();

    this->factors = THEKERNEL->config->value(leveling_strategy_checksum, lsq_delta_calibration_strategy_checksum, factors_checksum)->by_default(6)->as_number();
    
    // the initial height above the bed we stop the intial move down after home to find the bed
    // this should be a height that is enough that the probe will not hit the bed and is an offset from max_z (can be set to 0 if max_z takes into account the probe offset)
    this->initial_height= THEKERNEL->config->value(leveling_strategy_checksum, lsq_delta_calibration_strategy_checksum, initial_height_checksum)->by_default(10)->as_number();
    
    return true;
}

bool LSQDeltaCalibrationStrategy::handleGcode(Gcode *gcode)
{
    if (gcode->has_g) {
        // G code processing
        if (gcode->g == 32) { // LSQ auto calibration for delta
            // first wait for an empty queue i.e. no moves left
            THEKERNEL->conveyor->wait_for_idle();
            if (!calibrate(gcode)) {
                gcode->stream->printf("Calibration failed to complete, probe not triggered\r\n");
                return true;
            }
            gcode->stream->printf("Calibration complete, save settings with M500\r\n");
            return true;
        }
        if(gcode->g == 29){
            // first wait for an empty queue i.e. no moves left
            THEKERNEL->conveyor->wait_for_idle();
            if (!manual_probe(gcode)) {
                gcode->stream->printf("G29 failed to complete, probe not triggered\r\n");
                return true;
            }
            gcode->stream->printf("G29 complete\r\n");
            return true;
        }

    } else if(gcode->has_m) {
        // handle mcodes
    }

    return false;
}

bool LSQDeltaCalibrationStrategy::manual_probe(Gcode *gcode)
{
    int sample_count = this->sample_count;
    if (gcode->has_letter('S')) {
        sample_count = gcode->get_value('S');
    }
    
    if(sample_count > MAX_DELTA_PROBE_POINTS)
        sample_count = MAX_DELTA_PROBE_POINTS;  //extremely basic coerce to protect against out-of-bounds
    
    float probe_radius = this->probe_radius;
    if (gcode->has_letter('R')) {
        probe_radius = gcode->get_value('R');
    }
    
    float probe_heights[MAX_DELTA_PROBE_POINTS];
    probe_bed(sample_count, probe_radius, probe_heights, gcode);
    
    float sumOfSquares = 0;
    for (int i = 0; i < sample_count; i++) {
        sumOfSquares += (probe_heights[i] * probe_heights[i]);
    }
    gcode->stream->printf("RMS Error: %.4f\r\n", sqrtf(sumOfSquares/sample_count));
    
    return true;
}

/* Run the David Crocker least squares calibration algorithm for delta calibration.

   The David Crocker algorithm calculates the derivative, or sensitivity,
   of each component at each sample point.  This is used to calculate a 
   normal matrix for each component.  This normal matrix, combined with the
   height errors from the samples, is used to calculate a solution matrix
   for each factor that is tuned.
   
   This approach simplifies delta kinematics: individual tower radius 
   adjustments are mathematically identical to a new circle with a single 
   delta radius for all towers, with slight adjustments to each tower angle.
   Furthermore, one of the 3 tower angles is assumed to be '0', with all 
   other tower angles relative to the first.  This reduces an 11-variable
   equation (DR, RL, 3 endstops, 3 delta radius adjustments, 3 tower angle
   adjustments) to a mathematically identical, bunt sipler, equation with
   7 variables (DR, RL, 3 endstops, 2 tower angles).
   
   The usage is as follows:

   G32 (Fx) (Sx) (Rx.xx) (K)
        Fx : Number of factors.  Valid values are 3, 4, 6, 7.
             Default is 6.
             3 = calibrate endstops for towers A,B,C
             4 = calibrate endstops and delta radius
             6 = calibrate endstops, delta radius, and tower angles for A,B
             7 = calibrate endstops, delta radius, tower angles, and rod length
             NOTE: autocalibration of rod length generally should be avoided, 
			 if your object scale is uniformly wrong in X & Y, tweak rod length manually.
        Sx : Number of sample points. 
		     Should be a multiple of 6 plus 1.  
             Default is 7. Max is 13. Developers can try increasing 
			 MAX_DELTA_PROBE_POINTS and compiling to see if higher points work.
        Rx.xx : Maximum probe radius.  Default is 100.0.
             Probes always include the center.  
			 Probes are performed in a series of circles, 6 points per circle, 
			 varying the radius linearly from 0 to the probe radius.
     
 */ 
bool LSQDeltaCalibrationStrategy::calibrate(Gcode *gcode)
{
    int factors = this->factors;
    if (gcode->has_letter('F')) {
        factors = gcode->get_value('F');
        if (factors != 3 && factors != 4 && factors != 6 && factors != 7) {
            gcode->stream->printf("Number of factors for LSQ calibration is incorrect--must be 3, 4, 6, or 7\r\n");
            return false;
        }
    }
    
    int sample_count = this->sample_count;
    if (gcode->has_letter('S')) {
        sample_count = gcode->get_value('S');
    }
    
    if(sample_count > MAX_DELTA_PROBE_POINTS)
        sample_count = MAX_DELTA_PROBE_POINTS;  //extremely basic coerce to protect against out-of-bounds
    
    float probe_radius = this->probe_radius;
    if (gcode->has_letter('R')) {
        probe_radius = gcode->get_value('R');
    }
    
    bool keep = false;
    if (gcode->has_letter('K')) {
        keep = true;
    }
    
    gcode->stream->printf("Starting calibration. Factors: %i, sample count: %i, radius: %.4f, keep existing calibration? %i\r\n", factors, sample_count, probe_radius, keep);
    
    return calibrate(factors, sample_count, probe_radius, keep, gcode);
}

bool LSQDeltaCalibrationStrategy::calibrate(int num_factors, int sample_count, float probe_radius, bool keep, Gcode *gcode)
{
    float cartesian_mm[3];
    //std::array<float, 3> actuator_mm;
	ActuatorCoordinates actuator_mm;
    std::vector<std::array<float, 3>> all_actuator_mm(sample_count);
    BaseSolution::arm_options_t options;
    
    // Get or reset starting point for options, depending on 'keep'
    float trimx = 0.0F, trimy = 0.0F, trimz = 0.0F;
    if (!keep) {
        // zero trim values
        if (!set_trim(0, 0, 0)) return false;

        options['A'] = 0;
        options['B'] = 0;
        options['C'] = 0;
        options['D'] = 0;
        options['E'] = 0;
        options['F'] = 0;
        THEKERNEL->robot->arm_solution->set_optional(options);
    } else {
        // get current trim, and continue from that
        if (get_trim(trimx, trimy, trimz)) {
            gcode->stream->printf("Current Trim X: %.3f, Y: %.3f, Z: %.3f\r\n", trimx, trimy, trimz);
        } else {
            gcode->stream->printf("Could not get current trim, are endstops enabled?\r\n");
            return false;
        }
    }
    if (THEKERNEL->robot->arm_solution->get_optional(options, true)) {
        gcode->stream->printf("    Rod length: %.4f\r\n", options['L']);
        gcode->stream->printf("  Delta radius: %.4f\r\n", options['R']);
        gcode->stream->printf("Tower A radius: %.4f\r\n", options['R'] + options['A']);
        gcode->stream->printf("Tower A  angle: %.4f\r\n", 210 + options['D']);
        gcode->stream->printf("Tower B radius: %.4f\r\n", options['R'] + options['B']);
        gcode->stream->printf("Tower B  angle: %.4f\r\n", 330 + options['E']);
        gcode->stream->printf("Tower C radius: %.4f\r\n", options['R'] + options['C']);
        gcode->stream->printf("Tower C  angle: %.4f\r\n", 90 + options['F']);
    }

    // Sample bed Z height
    float probe_heights[MAX_DELTA_PROBE_POINTS];
    probe_bed(sample_count, probe_radius, probe_heights, gcode);
    
    // Remember actuator heights for probe points
    for (int i = 0; i < sample_count; i++) {
        get_probe_point(i, cartesian_mm, sample_count, probe_radius);
        cartesian_mm[2] = probe_heights[i];
        THEKERNEL->robot->arm_solution->cartesian_to_actuator(cartesian_mm, actuator_mm);
        all_actuator_mm[i][0] = actuator_mm[0];
        all_actuator_mm[i][1] = actuator_mm[1];
        all_actuator_mm[i][2] = actuator_mm[2];
        THEKERNEL->call_event(ON_IDLE);
    }

    // The amount of Z height, at each probe point, changed by the solution
    float corrections[MAX_DELTA_PROBE_POINTS];
    float initialSumOfSquares = 0;
    for (int i = 0; i < sample_count; i++) {
        corrections[i] = 0;
        initialSumOfSquares += (probe_heights[i] * probe_heights[i]);
        THEKERNEL->call_event(ON_IDLE);
    }
    
    gcode->stream->printf("Starting RMS Error: %.4f\r\n", sqrtf(initialSumOfSquares/sample_count));
    
    // The derivative of height change for each point wrt each factor
    FixedMatrix<float, MAX_DELTA_PROBE_POINTS, NumDeltaFactors> derivative_matrix;
    // The normal matrix for each factor
	FixedMatrix<float, NumDeltaFactors, NumDeltaFactors + 1> normal_matrix;
    
    // The error residuals after 
    //std::vector<float> residuals(sample_count, 0);
    //float residuals[MAX_DELTA_PROBE_POINTS];
    
    // The solution--the adjustments for each of the factors
    //std::vector<float> solution(num_factors, 0);
    
    for (int iteration = 0; iteration < 2; iteration++) {
        for (int i = 0; i < sample_count; i++) {
            get_probe_point(i, cartesian_mm, sample_count, probe_radius);
            cartesian_mm[2] = probe_heights[i];
            THEKERNEL->robot->arm_solution->cartesian_to_actuator(cartesian_mm, actuator_mm);   //Already have this information in all_actuator_mm... use it?
            for (int j = 0; j < num_factors; j++) {
                derivative_matrix(i,j) = compute_derivative(j, cartesian_mm, actuator_mm);
            }
            THEKERNEL->call_event(ON_IDLE);
        }
        
        print_matrix("Derivative Matrix", derivative_matrix, sample_count, num_factors, gcode);
        
        for (int i = 0; i < num_factors; i++)
        {
            for (int j = 0; j < num_factors; j++)
            {
                float temp = 0;
                for (int k = 0; k < sample_count; k++)
                {
                    temp += derivative_matrix(k,i) * derivative_matrix(k,j);
                }
                normal_matrix(i,j) = temp;
            }
            float temp = 0;
            for (int k = 0; k < sample_count; k++) {
                temp += derivative_matrix(k,i) * -(corrections[k] + probe_heights[k]);
            }
            normal_matrix(i,num_factors) = temp;
            THEKERNEL->call_event(ON_IDLE);
        }
        
        print_matrix("Normal matrix", normal_matrix, num_factors, num_factors +1, gcode);
        
        float solution[NumDeltaFactors];
        normal_matrix.GaussJordan(solution, num_factors);
        THEKERNEL->call_event(ON_IDLE);
        
        print_matrix("Solved matrix", normal_matrix, num_factors, num_factors+1, gcode);
        print_vector("Solution", solution, num_factors, gcode);
        
        // compute residuals
//        float sumOfSquares = 0;
//        float residuals[MAX_DELTA_PROBE_POINTS];
//        for (int i = 0; i < sample_count; i++) {
//            residuals[i] = probe_heights[i];
//            for (int j = 0; j < num_factors; j++) {
//                residuals[i] += solution[j] * derivative_matrix(i,j);
//                THEKERNEL->call_event(ON_IDLE);
//            }
//            sumOfSquares += residuals[i] * residuals[i];
//        }
//        stream->printf("Before RMS Error: %.3f\r\n", sqrtf(sumOfSquares));
        
        // apply solution
        trimx += -solution[0];
        trimy += -solution[1];
        trimz += -solution[2];
        if(!set_trim(trimx, trimy, trimz)) return false;
        if(num_factors >= 4){
            options['R'] += solution[3];
            if(num_factors >= 6){
                options['D'] += solution[4];
                options['E'] += solution[5];
                if(num_factors == 7){
                    options['L'] += solution[6];
                }
            }
        }
        
        THEKERNEL->robot->arm_solution->set_optional(options);
        gcode->stream->printf("X: %.4f, Y: %.4f, Z: %.4f, R: %.4f, D: %.4f, E: %.4f, F: 0.000, L: %.4f\r\n", trimx, trimy, trimz, options['R'], options['D'], options['E'], options['L']);
        THEKERNEL->call_event(ON_IDLE);
        
        
        // compute expected residuals
        float sumOfSquares = 0;
        float residuals[MAX_DELTA_PROBE_POINTS];
        for (int i = 0; i < sample_count; i++) {
            actuator_mm[0] = all_actuator_mm[i][0] + -trimx;
            actuator_mm[1] = all_actuator_mm[i][1] + -trimy;
            actuator_mm[2] = all_actuator_mm[i][2] + -trimz;
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            residuals[i] = cartesian_mm[2];
            corrections[i] = (cartesian_mm[2] - probe_heights[i]);
            sumOfSquares += residuals[i] * residuals[i];
            THEKERNEL->call_event(ON_IDLE);
        }
        gcode->stream->printf("RMS Error: %.4f\r\n", sqrtf(sumOfSquares/sample_count));
        
        print_vector("Corrections", corrections, sample_count, gcode);
        print_vector("Residuals", residuals, sample_count, gcode);
    }

    return true;
}

void LSQDeltaCalibrationStrategy::print_matrix(const char* s, const MathMatrix<float>& m, int maxRows, int maxCols, Gcode *gcode) {
    gcode->stream->printf("%s:\r\n", s);
    for (int i = 0; i < maxRows; i++) {
        for (int j = 0; j < maxCols; j++) {
            if (j == 0) {
                gcode->stream->printf("[");
            } else {
                gcode->stream->printf(",");
            }
            gcode->stream->printf("%.4f", m(i,j));
        }
        gcode->stream->printf("]\r\n");
    }
}

void LSQDeltaCalibrationStrategy::print_vector(const char* s, const float *v, int size, Gcode *gcode){
    gcode->stream->printf("%s:\r\n", s);
    for (int i = 0; i < size; i++) {
        if (i == 0) {
            gcode->stream->printf("[");
        } else {
            gcode->stream->printf(",");
        }
        gcode->stream->printf("%.4f", v[i]);
        THEKERNEL->call_event(ON_IDLE);
    }
    gcode->stream->printf("]\r\n");
}

bool LSQDeltaCalibrationStrategy::set_trim(float x, float y, float z)
{
    float t[3] {x, y, z};
    return PublicData::set_value( endstops_checksum, trim_checksum, t);
}

bool LSQDeltaCalibrationStrategy::get_trim(float &x, float &y, float &z)
{
    void *returned_data;
    bool ok = PublicData::get_value( endstops_checksum, trim_checksum, &returned_data );

    if (ok) {
        float *trim = static_cast<float *>(returned_data);
        x = trim[0];
        y = trim[1];
        z = trim[2];
        return true;
    }
    return false;
}

bool LSQDeltaCalibrationStrategy::probe_bed(int sample_count, float probe_radius, float* probe_heights, Gcode *gcode) {
    
	float mm;
    float bedht= findBed();
    if(isnan(bedht)) return false;
    gcode->stream->printf("Initial Bed height is %.4f mm\n", bedht);

    // move to start position
    zprobe->home();
    zprobe->coordinated_move(NAN, NAN, -bedht, zprobe->getFastFeedrate(), true); // do a relative move from home to the point above the bed

    float cartesian_mm[3];
    for (int i = 0; i < sample_count; i++) {
        get_probe_point(i, cartesian_mm, sample_count, probe_radius);
        if (!zprobe->doProbeAt(mm, cartesian_mm[0], cartesian_mm[1])) return false;
        probe_heights[i] = zprobe->getProbeHeight() - mm;
        gcode->stream->printf("X: %.4f, Y: %.4f, Z: %.4f\r\n", cartesian_mm[0], cartesian_mm[1], probe_heights[i]);
        THEKERNEL->call_event(ON_IDLE);
    }
    
    return true;
}

float LSQDeltaCalibrationStrategy::findBed()
{
    // home
    zprobe->home();
    
    // move to an initial position fast so as to not take all day, we move down max_z - initial_height, which is set in config, default 10mm
    float deltaz= zprobe->getMaxZ() - initial_height;
    zprobe->coordinated_move(NAN, NAN, -deltaz, zprobe->getFastFeedrate(), true);
    
    // find bed, run at slow rate so as to not hit bed hard
    float mm;
    if(!zprobe->run_probe(mm, false)) return NAN;
    
    return mm + deltaz - zprobe->getProbeHeight(); // distance to move from home to 5mm above bed
}

void LSQDeltaCalibrationStrategy::get_probe_point(int sample_number, float cartesian_mm[], int sample_count, float probe_radius) {
    // The last sample is always the center
    if (sample_number == (sample_count-1)) {
        cartesian_mm[0] = 0;
        cartesian_mm[1] = 0;
    } else {
        int circle_count = floor((sample_count - 1) / 6);
        // Work from the outside in...this is the circle #
        int circle = circle_count - floor(sample_number / 6);
        // Each successive circle gets closer to the center
        float radius = (probe_radius * circle) / floor((sample_count - 1) / 6);
        // Probe points are at 0, 60, 120, 180, 240, and 300 degress for odd numbered circles
        //              and at 30, 90, 150, 210, 270, and 330 degrees for even numbered circles
        float angle = PI * (sample_number % 6) / 3.0;
        if ((circle % 2) == 0 || sample_count == 7) {
            angle += (PI / 6.0);
        }
        
        cartesian_mm[0] = cosf(angle) * radius;
        cartesian_mm[1] = sinf(angle) * radius;
    }
}

float LSQDeltaCalibrationStrategy::compute_derivative(int factor, float cartesian_mm[], ActuatorCoordinates actuator_mm) {
    float perturb = 0.2;
    float zLo = 0;
    float zHi = 0;
    float original;
    BaseSolution::arm_options_t options;
    
    if (factor < 0) {
        // error
        return 0;
    } else if (factor < 3) {
        // perturb endstops
        original = actuator_mm[factor];
        // DSK: TODO: check that this works--it's opposite the DeltaSim
        actuator_mm[factor] = original + perturb; // endstop offsets have inverted effects
        THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
        zHi = cartesian_mm[2];
        
        actuator_mm[factor] = original - perturb; // endstop offsets have inverted effects
        THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
        zLo = cartesian_mm[2];
        
        actuator_mm[factor] = original;
    } else if (factor < 7) {
        // perturb option
        char option = '-';
        switch (factor) {
            case 3:
                option = 'R'; // Delta radius adjustment
                break;
            case 4:
                option = 'D'; // Angle adjustment for tower A
                break;
            case 5:
                option = 'E'; // Angle adjustment for tower B
                break;
            case 6:
                option = 'L'; // Rod length adjustment
                break;
        }
        if (THEKERNEL->robot->arm_solution->get_optional(options, true)) {
            original = options[option];
            
            options[option] = original + perturb;
            THEKERNEL->robot->arm_solution->set_optional(options);
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            zHi = cartesian_mm[2];
            
            options[option] = original - perturb;
            THEKERNEL->robot->arm_solution->set_optional(options);
            THEKERNEL->robot->arm_solution->actuator_to_cartesian(actuator_mm, cartesian_mm);
            zLo = cartesian_mm[2];

            options[option] = original;
            THEKERNEL->robot->arm_solution->set_optional(options);
        }
    } else {
        // error
        return 0;
    }
    return (zHi - zLo) / (2 * perturb);
}
