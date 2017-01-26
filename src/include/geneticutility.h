void calcPointInsideBoxInfo(const float * const pts, const float * const bxs,
                            const bool * const inout,
                            unsigned * res,
                            unsigned n, unsigned M, unsigned d);

float* readInputFile(const char* filename, unsigned& m, unsigned& d, unsigned& n);

