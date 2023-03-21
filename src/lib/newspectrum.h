#pragma once

#include "lib/spectrum.h"
#include "vec3.h"
#include <array>

// I went through section 5.1, Spectral Representation, from the book "Physically Based
// Rendering," and got a general idea of how to implement the new spectrum class.
//
// https://pbr-book.org/3ed-2018/Color_and_Radiometry/Spectral_Representation
//
// In order to render the prism-like effect, the spectrum cannot just have the RGB three numbers.
//
// We can see the rainbow effect because light at different wavelengths as different refractive
// indices when they travel through the glass.
//
// Although all the components in a white light enter the glass have the same incoming direction,
// because they have different refractive indexes, their outgoing direction are all different.
// This is how a prism can disperse light.
//
// Instead of storing RGB three values as a vector like in the old spectrum class, the new spectrum
// class needs to support more color channels. In our case, we have an array with 60 elements.
// Each element corresponds to the light intensity at a particular wavelength.
//
// The wavelength starts from 400 nm and ends at 700 nm with a interval of 5 nm.
static const float Lambda_Start = 400.0;
static const float Lambda_End = 700.0;
static const int nSpectrumSamples = 60;
static const float Lambda_Interval = (Lambda_End - Lambda_Start) / nSpectrumSamples;

// The book mentioned a method for converting RGBs to spectrum suggested by Smits(1999),
// which requires using some known smooth red, green, blue, cyan, magenta, and yellow spectrum.
// The length of those spectrum samples is 28.
static const int RGB2SpectrumSamples = 28;

// When converting a spectrum to XYZ, we need some known X(λ), Y(λ), and Z(λ) response curves.
// They are sampled at 5nm increments, so the array length is 60.
static const int nCIESamples = 60;

class NewSpectrum {
public:
    float c[nSpectrumSamples];
    // Initialization of NewSpectrum object with all 0 values.
    NewSpectrum() {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = 0;
    }
    // Initialization of NewSpectrum object with float a.
    NewSpectrum(float a) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = a;
    }
    // Initialization of NewSpectrum object with a known spectrum sample
    NewSpectrum(std::array<float, RGB2SpectrumSamples>& lambdas,
                std::array<float, RGB2SpectrumSamples>& values);
    // Add intensity at a specific wavelength to the current spectrum.
    void addValueAtLambda(float lambda, float c);
    // Return light intensity at a specific wavelength.
    float sampleAtLambda(float lambda);
    // Check if the spectrum is valid
    bool valid() {
        for(int i = 0; i < nSpectrumSamples; ++i) {
            if(std::isinf(c[i]) || std::isnan(c[i])) return false;
        }
        return true;
    }
    // Convert the current spectrum to an RGB value (Vec3)
    Vec3 Convert2RGB();
    // Operators
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
// inline Operators
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
//Convert an RGB value to a Spectrum
NewSpectrum Convert2Spectrum(float r, float g, float b);

//Convert an original Spectrum object to the new Spectrum object
NewSpectrum Old2NewSpectrum(Spectrum old_spectrum);
//Convert a new Spectrum object to the original Spectrum object
Spectrum New2OldSpectrum(NewSpectrum new_spectrum);
