#pragma once

#include "lib/spectrum.h"
#include "vec3.h"
#include <array>

static const float Lamda_Start = 380.0;
static const float Lamda_End = 780.0;
static const int nSpectrumSamples = 80;
static const float Lamda_Interval = (Lamda_End - Lamda_Start) / nSpectrumSamples;
static const int RGB2SpectrumSamples = 32;
static const int nCIESamples = 80;
class NewSpectrum {
public:
    float c[nSpectrumSamples];

    NewSpectrum() {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = 0;
    }

    NewSpectrum(float a) {
        for(int i = 0; i < nSpectrumSamples; ++i) c[i] = a;
    }

    NewSpectrum(float r, float g, float b);

    NewSpectrum(std::array<float, RGB2SpectrumSamples>& lambdas,
                std::array<float, RGB2SpectrumSamples>& values);

    void addValueAtLambda(float lambda, float c);
    float sampleAtLambda(float lambda);

    bool valid() {
        for(int i = 0; i < nSpectrumSamples; ++i) {
            if(std::isinf(c[i]) || std::isnan(c[i])) return false;
        }
        return true;
    }

    Vec3 Convert2RGB();

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

NewSpectrum Convert2Spectrum(float r, float g, float b);
NewSpectrum Old2NewSpectrum(Spectrum old_spectrum);
Spectrum New2OldSpectrum(NewSpectrum new_spectrum);
