#ifndef IMAGINE_FFTW_H
#define IMAGINE_FFTW_H

#include <array>
#include <cstddef>
#include <mutex>
#include <new>

#include <fftw3.h>

namespace imagine {

inline std::mutex &fftw_planner_mutex() {
    static std::mutex m;
    return m;
}

class FFTWWorkspace {
public:
    explicit FFTWWorkspace(const std::array<int, 3> &shape) : shape_(shape), padded_z_(2 * (shape[2] / 2 + 1)) {
        data_ = fftw_alloc_real(padded_size());
        if (!data_)
            throw std::bad_alloc();
        std::lock_guard<std::mutex> lock(fftw_planner_mutex());
        r2c_ = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], data_, complex(), FFTW_ESTIMATE);
        c2r_ = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], complex(), data_, FFTW_ESTIMATE);
    }

    ~FFTWWorkspace() {
        {
            std::lock_guard<std::mutex> lock(fftw_planner_mutex());
            fftw_destroy_plan(c2r_);
            fftw_destroy_plan(r2c_);
        }
        fftw_free(data_);
    }

    FFTWWorkspace(const FFTWWorkspace &) = delete;
    FFTWWorkspace &operator=(const FFTWWorkspace &) = delete;

    double *real() { return data_; }
    fftw_complex *complex() { return reinterpret_cast<fftw_complex *>(data_); }
    const std::array<int, 3> &shape() const { return shape_; }
    std::array<int, 3> padded_shape() const { return {shape_[0], shape_[1], padded_z_}; }
    std::size_t padded_size() const { return std::size_t(shape_[0]) * shape_[1] * padded_z_; }
    std::size_t size() const { return std::size_t(shape_[0]) * shape_[1] * shape_[2]; }

    void forward() { fftw_execute(r2c_); }
    void backward() { fftw_execute(c2r_); }

    void copy_unpadded(double *out) const {
        const std::size_t rows = std::size_t(shape_[0]) * shape_[1];
        for (std::size_t r = 0; r < rows; ++r)
            for (int k = 0; k < shape_[2]; ++k)
                out[r * shape_[2] + k] = data_[r * padded_z_ + k];
    }

private:
    std::array<int, 3> shape_;
    int padded_z_;
    double *data_ = nullptr;
    fftw_plan r2c_ = nullptr;
    fftw_plan c2r_ = nullptr;
};

}

#endif
