#ifndef _LSQDELTALEVELINGSTRATEGY
#define _LSQDELTALEVELINGSTRATEGY

#include <tuple>
#include <vector>
#include "LevelingStrategy.h"
#include "Matrix.h"
#include "ActuatorCoordinates.h"

#define lsq_delta_calibration_strategy_checksum CHECKSUM("lsq-delta-calibration")

class StreamOutput;

class LSQDeltaCalibrationStrategy : public LevelingStrategy
{
public:
    LSQDeltaCalibrationStrategy(ZProbe* zprobe) : LevelingStrategy(zprobe) {};
    ~LSQDeltaCalibrationStrategy() {};
    bool handleGcode(Gcode* gcode);
    bool handleConfig();
private:
    bool set_trim(float x, float y, float z);
    bool get_trim(float &x, float &y, float &z);
    float compute_derivative(int factor, float cartesian_mm[], ActuatorCoordinates actuator_mm);
    void get_probe_point(int sample_number, float cartesian_mm[], int sample_count, float probe_radius);
    bool probe_bed(int sample_count, float probe_radius, float* probe_heights, Gcode *gcode);
    float findBed();
    bool manual_probe(Gcode *gcode);
    bool calibrate(Gcode* gcode);
    bool calibrate(int numFactors, int sample_count, float probe_radius, bool keep, Gcode *gcode);
    void print_matrix(const char* s, const MathMatrix<float>& m, int maxRows, int maxCols, Gcode *gcode);
    void print_vector(const char* s, const float *v, int size, Gcode *gcode);
    
    float probe_radius;
    int sample_count;
    int factors;
    float initial_height;
};

#endif
