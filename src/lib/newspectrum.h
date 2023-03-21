#pragma once

#include "lib/spectrum.h"
#include "vec3.h"
#include <array>

// I went through the section 5.1 Spectral Representation from the book "Physically Based
// Rendering", and get a general idea on how to implement the new spectrum class.
//
// https://pbr-book.org/3ed-2018/Color_and_Radiometry/Spectral_Representation
//
// In order to render the prism-like effect, the spectrum cannot just have the RGB three numbers.
//
// The reason we can see the rainbow is because light at different wavelength have different
// refractive index when they travel through teh glass.
//
// Although all the components in a white light enter the glass have the same incoming direction,
// but because they have different refractive index, their out direction are all different. This is
// how a prism can disperse light.
//
// Instead storing RGB three values as a vector like in the old spectrum class, the new spectrum
// class needs to support more color channels. In our case, we have a array with 60 elements, each
// elements corresponds to the light intensity at a certain wavelength.
//
// The wavelength starts from 400 nm and ends at 700 nm with a interval of 5 nm.
static const float Lambda_Start = 400.0;
static const float Lambda_End = 700.0;
static const int nSpectrumSamples = 60;
static const float Lambda_Interval = (Lambda_End - Lambda_Start) / nSpectrumSamples;

// The book mentioned a method for converting RGBs to spectrum suggested by Smits(1999) which
// requires the use of some known smooth red, green, blue, cyan, magenta and yellow spectrum. The
// length of those spectrum sample are 28.
static const int RGB2SpectrumSamples = 28;

// When converting spectrum to XYZ, we need some known X(λ), Y(λ), Z(λ) response curves. They are
// sampled at 5nm increments, so the array length is 60.
static const int nCIESamples = 60;

class NewSpectrum {
public:
    float c[nSpectrumSamples];
    // initialization of NewSpectrum object with all 0 value.
    NewSpectrum() {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = 0;
    }
    // initialization of NewSpectrum object with float a.
    NewSpectrum(float a) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = a;
    }
    // initialization of NewSpectrum object with a known spectrum sample
    NewSpectrum(std::array<float, RGB2SpectrumSamples>& lambdas,
                std::array<float, RGB2SpectrumSamples>& values);
    // add intensity at certain wavelength to current spectrum
    void addValueAtLambda(float lambda, float c);
    // return light intensity at certain wavelength
    float sampleAtLambda(float lambda);
    // check if the spectrum is valid
    bool valid() {
        for(int i = 0; i < nSpectrumSamples; ++i) {
            if(std::isinf(c[i]) || std::isnan(c[i])) return false;
        }
        return true;
    }
    // convert current spectrum to RGB value (Vec3)
    Vec3 Convert2RGB();
    // operators
    NewSpectrum& operator+=(const NewSpectrum& s2) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] += s2.c[i];
        return *this;
    }
    NewSpectrum& operator+=(const float a) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] += a;
        return *this;
    }

    NewSpectrum& operator*=(const NewSpectrum& s2) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] *= s2.c[i];
        return *this;
    }
    NewSpectrum& operator*=(const float a) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] *= a;
        return *this;
    }
};
// operators
inline NewSpectrum operator-(const NewSpectrum& s1, const NewSpectrum& s2) {
    NewSpectrum ret;
    for(int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = s1.c[i] - s2.c[i];
    return ret;
}
inline NewSpectrum operator*(const NewSpectrum& s1, const NewSpectrum& s2) {
    NewSpectrum ret;
    for(int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = s1.c[i] * s2.c[i];
    return ret;
}
//convert RGB values to spectrum
NewSpectrum Convert2Spectrum(float r, float g, float b);

//convert original Spectrum object to new Spectrum object
NewSpectrum Old2NewSpectrum(Spectrum old_spectrum);
//convert new Spectrum object to original Spectrum object
Spectrum New2OldSpectrum(NewSpectrum new_spectrum);
