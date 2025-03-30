/*
 * DISTRHO Plugin Framework (DPF)
 * Copyright (C) 2012-2021 Filipe Coelho <falktx@falktx.com>
 * Copyright (C) 2019-2021 Jean Pierre Cimalando <jp-dev@inbox.ru>
 *
 * Permission to use, copy, modify, and/or distribute this software for any purpose with
 * or without fee is hereby granted, provided that the above copyright notice and this
 * permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD
 * TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN
 * NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
 * IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include "DistrhoUI.hpp"

#include "Artwork.hpp"
#include "DemoWidgetBanner.hpp"
#include "DemoWidgetClickable.hpp"
#include  "EditableText.hpp"
#include <cstdio>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

// Add this class to FildesUI.cpp before the FildesUI class
class TransferFunctionVisualizer {
public:
    TransferFunctionVisualizer() 
        : mNumerator(1, 0.0)
        , mDenominator(1, 1.0)
        , mResponseCache()
        , mFrequencies()
        , mMagnitudes()
        , mPhases()
        , mSampleRate(44100.0)
        , mNeedsUpdate(true)
        , mZeros()
        , mPoles()
    {
        // Initialize with default frequency points (logarithmic scale)
        updateFrequencyPoints();
    }

    void setNumerator(const std::vector<double>& coeffs) {
        mNumerator = coeffs;
        mNeedsUpdate = true;
        findZeros();

        // Debug output
        std::cout << "Numerator set to: ";
        for (double c : mNumerator) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }

    void setDenominator(const std::vector<double>& coeffs) {
        mDenominator = coeffs;
        mNeedsUpdate = true;
        findPoles();

        // Debug output
        std::cout << "Denominator set to: ";
        for (double c : mDenominator) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }

    void setSampleRate(double sampleRate) {
        mSampleRate = sampleRate;
        updateFrequencyPoints();
        mNeedsUpdate = true;
    }

    void updateFrequencyPoints() {
        // Create logarithmic frequency scale from 20Hz to Nyquist
        const double nyquist = mSampleRate / 2.0;
        const double minFreq = 20.0;
        const int numPoints = 200;
        
        mFrequencies.resize(numPoints);
        
        for (int i = 0; i < numPoints; ++i) {
            const double t = static_cast<double>(i) / (numPoints - 1);
            // Logarithmic scale from minFreq to nyquist
            mFrequencies[i] = minFreq * std::pow(nyquist / minFreq, t);
        }
        
        mMagnitudes.resize(numPoints);
        mPhases.resize(numPoints);
        mNeedsUpdate = true;
    }

    void calculateResponse() {
        if (!mNeedsUpdate) {
            return;
        }
        
        const size_t numFreqs = mFrequencies.size();
        
        for (size_t i = 0; i < numFreqs; ++i) {
            double frequency = mFrequencies[i];
            double omega = 2.0 * M_PI * frequency / mSampleRate;
            
            std::complex<double> numerator(0.0, 0.0);
            std::complex<double> denominator(0.0, 0.0);
            
            // Calculate numerator
            for (size_t n = 0; n < mNumerator.size(); ++n) {
                double angle = -omega * static_cast<double>(n);
                std::complex<double> z(cos(angle), sin(angle));
                numerator += mNumerator[n] * z;
            }
            
            // Calculate denominator
            for (size_t n = 0; n < mDenominator.size(); ++n) {
                double angle = -omega * static_cast<double>(n);
                std::complex<double> z(cos(angle), sin(angle));
                denominator += mDenominator[n] * z;
            }
            
            // Calculate transfer function
            std::complex<double> response = numerator / denominator;
            
            // Calculate magnitude in dB and phase
            mMagnitudes[i] = 20.0 * log10(std::abs(response) + 1e-6);
            mPhases[i] = std::arg(response);
        }
        
        mNeedsUpdate = false;
    }

    // Find roots of a polynomial using companion matrix method
    std::vector<std::complex<double>> findRoots(const std::vector<double>& coeffs) {
        std::vector<std::complex<double>> roots;
        
        // Handle special cases
        if (coeffs.empty()) return roots;
        
        // Remove leading zeros
        std::vector<double> c = coeffs;
        while (c.size() > 1 && std::abs(c.back()) < 1e-10) {
            c.pop_back();
        }
        
        // Special case: constant or linear polynomial
        if (c.size() <= 2) {
            if (c.size() == 2 && std::abs(c[1]) > 1e-10) {
                roots.push_back(std::complex<double>(-c[0] / c[1], 0.0));
            }
            return roots;
        }
        
        // Normalize the coefficients
        double lead = c.back();
        for (double& coef : c) {
            coef /= lead;
        }
        
        // Create companion matrix
        const int n = c.size() - 1;
        std::vector<std::vector<double>> companion(n, std::vector<double>(n, 0.0));
        
        // Fill the first row with negated coefficients
        for (int i = 0; i < n; ++i) {
            companion[0][i] = -c[i] / c[n];
        }
        
        // Fill the subdiagonal with ones
        for (int i = 1; i < n; ++i) {
            companion[i][i-1] = 1.0;
        }
        
        // TODO: Compute eigenvalues using a proper method
        // For now, return some placeholder roots for visualization
        
        // Generate some placeholder roots for demonstration
        // In a real implementation, compute the eigenvalues of the companion matrix
        
        // Calculate some reasonable roots based on polynomial degree
        double angle = 2.0 * M_PI / n;
        for (int i = 0; i < n; ++i) {
            double theta = i * angle;
            double radius = 0.8; // Inside unit circle for stable filters
            roots.push_back(std::complex<double>(radius * cos(theta), radius * sin(theta)));
        }
        
        return roots;
    }

    void findZeros() {
        // Find zeros of the transfer function (roots of numerator)
        mZeros = findRoots(mNumerator);
    }

    void findPoles() {
        // Find poles of the transfer function (roots of denominator)
        mPoles = findRoots(mDenominator);
    }

    // Update a pole or zero based on dragging
    void updatePoleZero(bool isPole, int index, const std::complex<double>& newValue) {
        std::vector<std::complex<double>>& roots = isPole ? mPoles : mZeros;
    
        if (index >= 0 && index < static_cast<int>(roots.size())) {
            roots[index] = newValue;
            
            // For complex roots, ensure conjugate pairs are maintained
            bool isComplex = std::abs(newValue.imag()) > 1e-10;
            if (isComplex) {
                // Find and update the conjugate pair if it exists
                for (size_t i = 0; i < roots.size(); ++i) {
                    if (i != index && std::abs(roots[i].real() - newValue.real()) < 1e-10 && 
                        std::abs(roots[i].imag() + newValue.imag()) < 1e-10) {
                        // Found the conjugate pair, update it
                        roots[i] = std::complex<double>(newValue.real(), -newValue.imag());
                        break;
                    }
                }
            }
            
            // Immediately update the corresponding coefficient array
            if (isPole) {
                std::vector<double> newDenominator = rootsToCoeffs(mPoles);
                
                // Normalize denominator so that z^0 term is 1.0
                if (!newDenominator.empty() && std::abs(newDenominator[0]) > 1e-10) {
                    double z0Coeff = newDenominator[0];
                    for (double& coeff : newDenominator) {
                        coeff /= z0Coeff;
                    }
                }
                
                mDenominator = newDenominator;
            } else {
                mNumerator = rootsToCoeffs(mZeros);
            }
            
            // Mark for response update
            mNeedsUpdate = true;
        }
    }

    // Convert roots back to polynomial coefficients
    std::vector<double> rootsToCoeffs(const std::vector<std::complex<double>>& roots) {
        // Special case: no roots
        if (roots.empty()) {
            return {1.0};
        }
        
        // Start with the term (z - r_0)
        std::vector<std::complex<double>> coeffs = {-roots[0], 1.0};
        
        // Multiply by each (z - r_i) term
        for (size_t i = 1; i < roots.size(); ++i) {
            std::vector<std::complex<double>> newCoeffs(coeffs.size() + 1, std::complex<double>(0.0, 0.0));
            
            // Multiply coeffs by (z - r_i)
            for (size_t j = 0; j < coeffs.size(); ++j) {
                // Term for z^j * z (increases the power)
                newCoeffs[j + 1] += coeffs[j];
                
                // Term for z^j * (-r_i) (same power)
                newCoeffs[j] -= coeffs[j] * roots[i];
            }
            
            coeffs = newCoeffs;
        }
        
        // Convert complex coefficients to real
        std::vector<double> realCoeffs(coeffs.size(), 0.0);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            // Handle very small values that might be numerical artifacts
            if (std::abs(coeffs[i].real()) < 1e-12) {
                realCoeffs[i] = 0.0;
            } else {
                realCoeffs[i] = coeffs[i].real();
            }
            
            // Print warning for non-negligible imaginary parts
            if (std::abs(coeffs[i].imag()) > 1e-8) {
                std::cout << "Warning: Significant imaginary component in coefficient: " 
                        << coeffs[i].imag() << std::endl;
            }
        }
        
        return realCoeffs;
    }

    // Update coefficients from poles and zeros
    bool updateCoefficientsFromRoots() {
        // Convert zeros to numerator coefficients
        std::vector<double> newNum = rootsToCoeffs(mZeros);
        
        // Convert poles to denominator coefficients
        std::vector<double> newDen = rootsToCoeffs(mPoles);
        
        // Normalize denominator so that z^0 term is 1.0
        if (!newDen.empty() && std::abs(newDen[0]) > 1e-10) {
            double z0Coeff = newDen[0];
            for (double& coeff : newDen) {
                coeff /= z0Coeff;
            }
        }
        
        // Check if coefficients have changed
        bool changed = false;
        
        if (newNum.size() != mNumerator.size()) {
            mNumerator = newNum;
            changed = true;
        } else {
            for (size_t i = 0; i < newNum.size(); ++i) {
                if (std::abs(newNum[i] - mNumerator[i]) > 1e-10) {
                    mNumerator = newNum;
                    changed = true;
                    break;
                }
            }
        }
        
        if (newDen.size() != mDenominator.size()) {
            mDenominator = newDen;
            changed = true;
        } else {
            for (size_t i = 0; i < newDen.size(); ++i) {
                if (std::abs(newDen[i] - mDenominator[i]) > 1e-10) {
                    mDenominator = newDen;
                    changed = true;
                    break;
                }
            }
        }
        
        if (changed) {
            mNeedsUpdate = true;
        }
        
        return changed;
    }

    const std::vector<double>& getFrequencies() const {
        return mFrequencies;
    }

    const std::vector<double>& getMagnitudes() const {
        return mMagnitudes;
    }

    const std::vector<double>& getPhases() const {
        return mPhases;
    }

    const std::vector<std::complex<double>>& getZeros() const {
        return mZeros;
    }

    const std::vector<std::complex<double>>& getPoles() const {
        return mPoles;
    }
    
    const std::vector<double>& getNumerator() const {
        return mNumerator;
    }
    
    const std::vector<double>& getDenominator() const {
        return mDenominator;
    }

    double getMinMagnitude() const {
        return *std::min_element(mMagnitudes.begin(), mMagnitudes.end());
    }

    double getMaxMagnitude() const {
        return *std::max_element(mMagnitudes.begin(), mMagnitudes.end());
    }

private:
    std::vector<double> mNumerator;
    std::vector<double> mDenominator;
    std::vector<std::complex<double>> mResponseCache;
    std::vector<double> mFrequencies;
    std::vector<double> mMagnitudes;
    std::vector<double> mPhases;
    std::vector<std::complex<double>> mZeros;
    std::vector<std::complex<double>> mPoles; 
    double mSampleRate;
    bool mNeedsUpdate;
};

START_NAMESPACE_DISTRHO

// We need a few classes from DGL.
using DGL_NAMESPACE::CairoGraphicsContext;
using DGL_NAMESPACE::CairoImage;
using DGL_NAMESPACE::CairoImageKnob;
using DGL_NAMESPACE::CairoImageSwitch;

// And from ourselves
using DGL_NAMESPACE::DemoWidgetBanner;
using DGL_NAMESPACE::EditableText;

class FildesUI : public UI,
                       public CairoImageKnob::Callback,
                       public CairoImageSwitch::Callback,
                       public DemoWidgetClickable::Callback,
                       public EditableText::Callback
{
    ScopedPointer<CairoImageKnob> fKnob;
    ScopedPointer<CairoImageSwitch> fButton;
    ScopedPointer<DemoWidgetClickable> fWidgetClickable;
    ScopedPointer<EditableText> fCoeffT, fCoeffB;
    ScopedPointer<DGL::SubWidget> fContainer;

    // Add the visualizer
    TransferFunctionVisualizer fVisualizer;
    
    // Store the current coefficients
    std::vector<double> fNumerator;
    std::vector<double> fDenominator;
    
    // Area for the frequency response plot
    DGL::Rectangle<int> fResponseArea;

    // Area for the pole-zero plot
    DGL::Rectangle<int> fPoleZeroArea;
    
    // Dragging state
    bool fIsDragging;
    bool fDraggingPole;  // true for pole, false for zero
    int fDraggedIndex;
    double fDragStartX, fDragStartY;

public:
    FildesUI()
    {
        freopen("output.txt","w",stdout);
        std::cout << "start\n" << std::flush;
        
        // Initialize with default filter coefficients
        parseCoefficients("0,0,0,1", fNumerator);
        parseCoefficients("1,0,0,0", fDenominator);
        fVisualizer.setNumerator(fNumerator);
        fVisualizer.setDenominator(fDenominator);
        
        CairoImage knobSkin;
        knobSkin.loadFromPNG(Artwork::knobData, Artwork::knobDataSize);

        fWidgetClickable = new DemoWidgetClickable(this);
        fWidgetClickable->setAbsolutePos(100, 100);
        fWidgetClickable->setSize(50, 50);
        fWidgetClickable->setCallback(this);
        fWidgetClickable->setId(kParameterTriState);
        fWidgetClickable->setVisible(false);

        fKnob = new CairoImageKnob(this, knobSkin);
        fKnob->setAbsolutePos(10, 100);
        fKnob->setSize(80, 80);
        fKnob->setCallback(this);
        fKnob->setId(kParameterKnob);
        fKnob->setVisible(false);

        CairoImage buttonOn, buttonOff;
        buttonOn.loadFromPNG(Artwork::buttonOnData, Artwork::buttonOnDataSize);
        buttonOff.loadFromPNG(Artwork::buttonOffData, Artwork::buttonOffDataSize);

        fButton = new CairoImageSwitch(this, buttonOff, buttonOn);
        fButton->setAbsolutePos(100, 160);
        fButton->setSize(60, 35);
        fButton->setCallback(this);
        fButton->setId(kParameterButton);
        fButton->setVisible(false);

        // Initialize in constructor
        fCoeffT = new EditableText(this, true);
        fCoeffT->setAbsolutePos(150, 50);
        fCoeffT->setSize(500, 30);
        fCoeffT->setCallback(this);

        fCoeffB = new EditableText(this, false);
        fCoeffB->setAbsolutePos(150, 100);
        fCoeffB->setSize(500, 30);
        fCoeffB->setCallback(this);
        
        // Set up the response and pole-zero areas
        fResponseArea = DGL::Rectangle<int>(50, 150, 700, 200);
        fPoleZeroArea = DGL::Rectangle<int>(500, 370, 250, 250);
        
        // Initialize dragging state
        fIsDragging = false;
        fDraggingPole = false;
        fDraggedIndex = -1;

        // Increase window size to accommodate both visualizations
        const double scaleFactor = getScaleFactor();
        setSize(800, 650);
        if (scaleFactor != 1.0)
            setSize(800 * scaleFactor, 650 * scaleFactor);
    }

protected:
    void onCairoDisplay(const CairoGraphicsContext& context)
    {
        cairo_t* const cr = context.handle;
        cairo_set_source_rgb(cr, 0.541, 0.576, 1);
        cairo_paint(cr);

        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); // Set text color
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 20);

        cairo_move_to(cr, 10, 30);
        cairo_show_text(cr, "Transfer Function:");
        
        cairo_move_to(cr, 10, 75);
        cairo_show_text(cr, "Numerator:");
        
        cairo_move_to(cr, 10, 125);
        cairo_show_text(cr, "Denominator:");

        cairo_move_to(cr, 10, 30);
        cairo_show_text(cr, "Transfer Function:"); // Static text

        drawFrequencyResponse(cr);

        drawPoleZeroPlot(cr);
    }

    void drawFrequencyResponse(cairo_t* cr)
    {
        // Calculate the response if needed
        fVisualizer.calculateResponse();
        
        // Prepare the plot area
        const int x = fResponseArea.getX();
        const int y = fResponseArea.getY();
        const int width = fResponseArea.getWidth();
        const int height = fResponseArea.getHeight();
        
        // Draw border and background
        cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
        cairo_rectangle(cr, x, y, width, height);
        cairo_fill_preserve(cr);
        cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
        cairo_set_line_width(cr, 1.0);
        cairo_stroke(cr);
        
        // Draw grid lines
        drawGrid(cr, x, y, width, height);
        
        // Get frequency response data
        const std::vector<double>& frequencies = fVisualizer.getFrequencies();
        const std::vector<double>& magnitudes = fVisualizer.getMagnitudes();
        
        if (frequencies.empty() || magnitudes.empty()) {
            return;
        }
        
        // Define magnitude range in dB (adjust as needed)
        const double minDb = -60.0;
        const double maxDb = 20.0;
        
        // Draw the magnitude response
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.8); // Blue for magnitude
        cairo_set_line_width(cr, 2.0);
        
        cairo_new_path(cr);
        
        for (size_t i = 0; i < frequencies.size(); ++i) {
            // Map frequency to x (logarithmic)
            const double logMinFreq = log10(20.0);
            const double logMaxFreq = log10(22050.0);
            const double logFreq = log10(frequencies[i]);
            const double normalizedX = (logFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            const double plotX = x + normalizedX * width;
            
            // Map magnitude to y (linear in dB)
            const double normalizedY = (magnitudes[i] - minDb) / (maxDb - minDb);
            const double clampedY = std::max(0.0, std::min(1.0, normalizedY));
            const double plotY = y + height - (clampedY * height);
            
            if (i == 0) {
                cairo_move_to(cr, plotX, plotY);
            } else {
                cairo_line_to(cr, plotX, plotY);
            }
        }
        
        cairo_stroke(cr);
        
        // Draw labels
        drawFrequencyLabels(cr, x, y, width, height);
        drawMagnitudeLabels(cr, x, y, width, height, minDb, maxDb);
    }

    void drawPoleZeroPlot(cairo_t* cr)
    {
        // Prepare the plot area
        const int x = fPoleZeroArea.getX();
        const int y = fPoleZeroArea.getY();
        const int size = fPoleZeroArea.getHeight(); // Square plot
        const int radius = size / 2 - 10;
        const int centerX = x + size / 2;
        const int centerY = y + size / 2;
        
        // Draw border and background
        cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
        cairo_rectangle(cr, x, y, size, size);
        cairo_fill_preserve(cr);
        cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
        cairo_set_line_width(cr, 1.0);
        cairo_stroke(cr);
        
        // Draw unit circle
        cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
        cairo_set_line_width(cr, 1.0);
        cairo_arc(cr, centerX, centerY, radius, 0, 2 * M_PI);
        cairo_stroke(cr);
        
        // Draw axes
        cairo_move_to(cr, centerX - radius - 5, centerY);
        cairo_line_to(cr, centerX + radius + 5, centerY);
        cairo_move_to(cr, centerX, centerY - radius - 5);
        cairo_line_to(cr, centerX, centerY + radius + 5);
        cairo_stroke(cr);
        
        // Draw poles (x)
        cairo_set_source_rgb(cr, 0.8, 0.0, 0.0); // Red for poles
        const std::vector<std::complex<double>>& poles = fVisualizer.getPoles();
        
        for (const auto& pole : poles) {
            // Map complex coordinate to plot coordinate
            const double plotX = centerX + pole.real() * radius;
            const double plotY = centerY - pole.imag() * radius;
            
            // Draw X
            const double crossSize = 6.0;
            cairo_set_line_width(cr, 2.0);
            cairo_move_to(cr, plotX - crossSize, plotY - crossSize);
            cairo_line_to(cr, plotX + crossSize, plotY + crossSize);
            cairo_move_to(cr, plotX - crossSize, plotY + crossSize);
            cairo_line_to(cr, plotX + crossSize, plotY - crossSize);
            cairo_stroke(cr);
        }
        
        // Draw zeros (o)
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.8); // Blue for zeros
        const std::vector<std::complex<double>>& zeros = fVisualizer.getZeros();
        
        for (const auto& zero : zeros) {
            // Map complex coordinate to plot coordinate
            const double plotX = centerX + zero.real() * radius;
            const double plotY = centerY - zero.imag() * radius;
            
            // Draw O
            const double circleSize = 6.0;
            cairo_set_line_width(cr, 2.0);
            cairo_arc(cr, plotX, plotY, circleSize, 0, 2 * M_PI);
            cairo_stroke(cr);
        }
        
        // Draw title
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 14);
        cairo_move_to(cr, x + 10, y + 20);
        cairo_show_text(cr, "Pole-Zero Plot");
        
        // Draw drag instruction
        cairo_set_font_size(cr, 10);
        cairo_move_to(cr, x + 10, y + size - 10);
        cairo_show_text(cr, "Drag poles (x) and zeros (o) to modify filter");
    }

    void drawGrid(cairo_t* cr, int x, int y, int width, int height)
    {
        // Draw horizontal grid lines (magnitude)
        cairo_set_source_rgba(cr, 0.5, 0.5, 0.5, 0.5);
        cairo_set_line_width(cr, 0.5);
        
        const int numMagLines = 9; // -60, -50, -40, -30, -20, -10, 0, 10, 20 dB
        for (int i = 0; i < numMagLines; ++i) {
            const double yPos = y + i * (height / (numMagLines - 1));
            cairo_move_to(cr, x, yPos);
            cairo_line_to(cr, x + width, yPos);
        }
        
        // Draw vertical grid lines (frequency)
        const double freqs[] = {20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000};
        const int numFreqs = sizeof(freqs) / sizeof(freqs[0]);
        
        for (int i = 0; i < numFreqs; ++i) {
            const double logMinFreq = log10(20.0);
            const double logMaxFreq = log10(22050.0);
            const double logFreq = log10(freqs[i]);
            const double normalizedX = (logFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            const double xPos = x + normalizedX * width;
            
            cairo_move_to(cr, xPos, y);
            cairo_line_to(cr, xPos, y + height);
        }
        
        cairo_stroke(cr);
    }

    void drawFrequencyLabels(cairo_t* cr, int x, int y, int width, int height)
    {
        // Draw frequency labels
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10);
        
        const double freqs[] = {20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000};
        const int numFreqs = sizeof(freqs) / sizeof(freqs[0]);
        
        for (int i = 0; i < numFreqs; ++i) {
            const double logMinFreq = log10(20.0);
            const double logMaxFreq = log10(22050.0);
            const double logFreq = log10(freqs[i]);
            const double normalizedX = (logFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            const double xPos = x + normalizedX * width;
            
            // Format frequency label
            std::string label;
            if (freqs[i] >= 1000) {
                label = std::to_string(static_cast<int>(freqs[i] / 1000)) + "k";
            } else {
                label = std::to_string(static_cast<int>(freqs[i]));
            }
            
            cairo_text_extents_t extents;
            cairo_text_extents(cr, label.c_str(), &extents);
            
            cairo_move_to(cr, xPos - extents.width / 2, y + height + 15);
            cairo_show_text(cr, label.c_str());
        }
        
        // X-axis label
        cairo_set_font_size(cr, 12);
        cairo_text_extents_t extents;
        const char* xLabel = "Frequency (Hz)";
        cairo_text_extents(cr, xLabel, &extents);
        cairo_move_to(cr, x + (width - extents.width) / 2, y + height + 30);
        cairo_show_text(cr, xLabel);
    }
    
    void drawMagnitudeLabels(cairo_t* cr, int x, int y, int width, int height, double minDb, double maxDb)
    {
        // Draw magnitude labels
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 10);
        
        const int numLabels = 9; // -60, -50, -40, -30, -20, -10, 0, 10, 20 dB
        for (int i = 0; i < numLabels; ++i) {
            const double db = minDb + i * ((maxDb - minDb) / (numLabels - 1));
            const double yPos = y + height - i * (height / (numLabels - 1));
            
            std::string label = std::to_string(static_cast<int>(db)) + " dB";
            
            cairo_text_extents_t extents;
            cairo_text_extents(cr, label.c_str(), &extents);
            
            cairo_move_to(cr, x - extents.width - 5, yPos + extents.height / 2);
            cairo_show_text(cr, label.c_str());
        }
        
        // Y-axis label
        cairo_set_font_size(cr, 12);
        cairo_save(cr);
        cairo_text_extents_t extents;
        const char* yLabel = "Magnitude (dB)";
        cairo_text_extents(cr, yLabel, &extents);
        cairo_move_to(cr, x - 40, y + (height + extents.width) / 2);
        cairo_rotate(cr, -M_PI / 2);
        cairo_show_text(cr, yLabel);
        cairo_restore(cr);
    }
    
    void parseCoefficients(const char* str, std::vector<double>& coeffs, bool isDenominator = false) {
        coeffs.clear();
    
        // Create a copy of the input string to tokenize
        char inputCopy[1024];
        strncpy(inputCopy, str, sizeof(inputCopy) - 1);
        inputCopy[sizeof(inputCopy) - 1] = '\0';
        
        // Parse comma-separated values into a temporary vector
        std::vector<double> tempCoeffs;
        char* token = strtok(inputCopy, ",");
        while (token != NULL) {
            double value = (token[0] == '\0') ? 0.0 : std::strtod(token, NULL);
            tempCoeffs.push_back(value);
            token = strtok(NULL, ",");
        }
        
        if (isDenominator) {
            // For denominator, add z^0 term with coefficient 1.0
            coeffs.push_back(1.0);
            // Then add the parsed coefficients (which are z^-1 and higher terms)
            coeffs.insert(coeffs.end(), tempCoeffs.begin(), tempCoeffs.end());
        } else {
            // For numerator, use the coefficients as is
            coeffs = tempCoeffs;
        }
        
        // Debug: print the parsed coefficients
        std::cout << "Parsed coefficients: ";
        for (double c : coeffs) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }

    // Format coefficients to string for display
    std::string formatCoefficients(const std::vector<double>& coeffs, bool isDenominator = false) {
        std::string result;
    
        // For denominator, we need to normalize to make z^0 term = 1 and skip it in the output
        std::vector<double> normalizedCoeffs;
        
        if (isDenominator && coeffs.size() > 0) {
            // Get the z^0 coefficient (should be the first element)
            double z0Coeff = coeffs[0];
            
            // Skip the z^0 term and normalize the rest by dividing by z0Coeff
            if (std::abs(z0Coeff) > 1e-10) {
                normalizedCoeffs.reserve(coeffs.size() - 1);
                for (size_t i = 1; i < coeffs.size(); ++i) {
                    normalizedCoeffs.push_back(coeffs[i] / z0Coeff);
                }
            } else {
                // If z^0 coefficient is close to zero, just use remaining coefficients as is
                normalizedCoeffs.assign(coeffs.begin() + 1, coeffs.end());
            }
        } else {
            // For numerator, use coefficients as is
            normalizedCoeffs = coeffs;
        }
        
        // Format the coefficients
        for (size_t i = 0; i < normalizedCoeffs.size(); ++i) {
            if (i > 0) {
                result += ",";
            }
            
            // Format with fixed precision - adjust the precision as needed
            char buffer[32];
            
            // Use fewer decimal places for values close to integers
            if (std::abs(std::round(normalizedCoeffs[i]) - normalizedCoeffs[i]) < 1e-10) {
                snprintf(buffer, sizeof(buffer), "%.0f", std::round(normalizedCoeffs[i]));
            } else {
                snprintf(buffer, sizeof(buffer), "%.6g", normalizedCoeffs[i]);
            }
            
            result += buffer;
        }
        
        return result;
    }

    void updateTransferFunction() {
        // Update the visualizer with current coefficients
        fVisualizer.setNumerator(fNumerator);
        fVisualizer.setDenominator(fDenominator);
        
        // Force a repaint to update the graph
        repaint();
    }

    void updateVisualizationFromCoefficients() {
        // Ensure the visualizer has the latest coefficients
        fVisualizer.setNumerator(fNumerator);
        fVisualizer.setDenominator(fDenominator);
        
        // Make sure poles and zeros are calculated
        fVisualizer.findPoles();
        fVisualizer.findZeros();
        
        // Force a repaint
        repaint();
    }

    // Check if a point is inside the pole-zero plot
    bool isPoleZeroAreaClick(int x, int y) {
        return fPoleZeroArea.contains(x, y);
    }
    
    // Find the closest pole or zero to a point
    bool findClosestPoleZero(int x, int y, bool& isPole, int& index) {
        const int centerX = fPoleZeroArea.getX() + fPoleZeroArea.getWidth() / 2;
        const int centerY = fPoleZeroArea.getY() + fPoleZeroArea.getHeight() / 2;
        const int radius = fPoleZeroArea.getHeight() / 2 - 10;
        
        // Convert screen coordinates to complex plane coordinates
        double complexX = static_cast<double>(x - centerX) / radius;
        double complexY = -static_cast<double>(y - centerY) / radius;
        
        // Get poles and zeros
        const std::vector<std::complex<double>>& poles = fVisualizer.getPoles();
        const std::vector<std::complex<double>>& zeros = fVisualizer.getZeros();
        
        double minDistance = 0.15; // Minimum distance to consider a hit (in complex plane units)
        double closest = 1000.0;
        
        // Check poles
        for (size_t i = 0; i < poles.size(); ++i) {
            double dx = poles[i].real() - complexX;
            double dy = poles[i].imag() - complexY;
            double distance = std::sqrt(dx*dx + dy*dy);
            
            if (distance < minDistance && distance < closest) {
                closest = distance;
                isPole = true;
                index = i;
            }
        }
        
        // Check zeros
        for (size_t i = 0; i < zeros.size(); ++i) {
            double dx = zeros[i].real() - complexX;
            double dy = zeros[i].imag() - complexY;
            double distance = std::sqrt(dx*dx + dy*dy);
            
            if (distance < minDistance && distance < closest) {
                closest = distance;
                isPole = false;
                index = i;
            }
        }
        
        return closest < minDistance;
    }
    
    // Convert screen coordinates to complex plane coordinates
    std::complex<double> screenToComplex(int x, int y) {
        const int centerX = fPoleZeroArea.getX() + fPoleZeroArea.getWidth() / 2;
        const int centerY = fPoleZeroArea.getY() + fPoleZeroArea.getHeight() / 2;
        const int radius = fPoleZeroArea.getHeight() / 2 - 10;
        
        double real = static_cast<double>(x - centerX) / radius;
        double imag = -static_cast<double>(y - centerY) / radius;
        
        // Limit to unit circle (with small margin)
        double magnitude = std::sqrt(real*real + imag*imag);
        if (magnitude > 0.95) {
            real *= 0.95 / magnitude;
            imag *= 0.95 / magnitude;
        }
        
        return std::complex<double>(real, imag);
    }
    
    // Update poles and zeros from dragging
    void updateCoefficientsFromDrag(int x, int y) {
        if (!fIsDragging || fDraggedIndex < 0) {
            return;
        }
        
        // Convert screen coordinates to complex plane
        std::complex<double> newValue = screenToComplex(x, y);
        
        // Update pole or zero
        fVisualizer.updatePoleZero(fDraggingPole, fDraggedIndex, newValue);
        
        // Get updated coefficients 
        std::vector<double> num = fVisualizer.getNumerator();
        std::vector<double> den = fVisualizer.getDenominator();
        
        // Update stored coefficients
        fNumerator = num;
        fDenominator = den;
        
        // Format coefficients for display
        std::string numStr = formatCoefficients(num, false);
        std::string denStr = formatCoefficients(den, true);
        
        // Debugging output for coefficient changes
        std::cout << "Updating text fields - Num: " << numStr << " Den: " << denStr << std::endl;
        
        // Directly update the text of EditableText widgets
        if (fCoeffT != nullptr) {
            // Force text update for numerator
            fCoeffT->setText(numStr);
        }
        
        if (fCoeffB != nullptr) {
            // Force text update for denominator
            fCoeffB->setText(denStr);
        }
        
        // Force repaint to ensure UI updates
        repaint();
        
        // Update the plugin state
        setState("topTF", numStr.c_str());
        setState("bottomTF", denStr.c_str());
    }

    // Mouse event handlers for dragging poles/zeros
    bool onMouse(const MouseEvent& ev) override {
        if (isPoleZeroAreaClick(ev.pos.getX(), ev.pos.getY())) {
            if (ev.button == 1) { // Left button
                if (ev.press) {
                    // Mouse down
                    bool isPole = false;
                    int index = -1;
                    
                    if (findClosestPoleZero(ev.pos.getX(), ev.pos.getY(), isPole, index)) {
                        // Start dragging
                        fIsDragging = true;
                        fDraggingPole = isPole;
                        fDraggedIndex = index;
                        fDragStartX = ev.pos.getX();
                        fDragStartY = ev.pos.getY();
                        return true;
                    }
                } else {
                    // Mouse up
                    if (fIsDragging) {
                        fIsDragging = false;
                        fDraggedIndex = -1;
                        return true;
                    }
                }
            }
        }
        
        return UI::onMouse(ev);
    }
    
    bool onMotion(const MotionEvent& ev) override {
        if (fIsDragging && fDraggedIndex >= 0) {
            // Handle dragging
            updateCoefficientsFromDrag(ev.pos.getX(), ev.pos.getY());
            return true;
        }
        
        return UI::onMotion(ev);
    }

    // we can use this if/when our resources are scalable, for now they are PNGs
    void onResize(const ResizeEvent& ev) override
    {
        UI::onResize(ev);

        const double scaleFactor = getScaleFactor();

        fWidgetClickable->setSize(50*scaleFactor, 50*scaleFactor);
        fWidgetClickable->setAbsolutePos(100*scaleFactor, 100*scaleFactor);

        fKnob->setSize(80*scaleFactor, 80*scaleFactor);
        fKnob->setAbsolutePos(10*scaleFactor, 100*scaleFactor);

        fButton->setSize(60*scaleFactor, 35*scaleFactor);
        fButton->setAbsolutePos(100*scaleFactor, 160*scaleFactor);
    
        // Update the response area
        fResponseArea = DGL::Rectangle<int>(
            50 * scaleFactor, 
            150 * scaleFactor, 
            700 * scaleFactor, 
            200 * scaleFactor
        );
    }

    void parameterChanged(const uint32_t index, const float value) override
    {
        switch (index)
        {
        case kParameterKnob:
            fKnob->setValue(value);
            break;
        case kParameterTriState:
            fWidgetClickable->setColorId(static_cast<int>(value + 0.5f));
            break;
        case kParameterButton:
            fButton->setDown(value > 0.5f);
            break;
        }
    }
    
    void focus(EditableText* current) override {
        if (current == fCoeffT) {
            fCoeffB->setFocused(false);
        } else if (current == fCoeffB) {
            fCoeffT->setFocused(false);
        } else {
            fCoeffB->setFocused(false);
            fCoeffT->setFocused(false);
        }
        repaint();
    }

    void setTF(std::string newText) override {
        std::cout << "setTF\n" << std::flush;
        if (newText[0] == 'T') {
            const char* coeffStr = &newText[1];
            parseCoefficients(coeffStr, fNumerator, false);
            updateTransferFunction();
            setState("topTF", &newText[1]);
        }
        else if (newText[0] == 'B') {
            const char* coeffStr = &newText[1];
            parseCoefficients(coeffStr, fDenominator, true);
            updateTransferFunction();
            setState("bottomTF", &newText[1]);
        }
    }

    void stateChanged(const char*, const char*) {
        return;
    }

    void demoWidgetClicked(DemoWidgetClickable*, const uint8_t colorId) override
    {
        setParameterValue(kParameterTriState, colorId);
    }

    void imageKnobDragStarted(CairoImageKnob*) override
    {
        editParameter(kParameterKnob, true);
    }

    void imageKnobDragFinished(CairoImageKnob*) override
    {
        editParameter(kParameterKnob, false);
    }

    void imageKnobValueChanged(CairoImageKnob*, const float value) override
    {
        setParameterValue(kParameterKnob, value);
    }

    void imageSwitchClicked(CairoImageSwitch*, bool down) override
    {
        setParameterValue(kParameterButton, down ? 1.f : 0.f);
    }
};

UI* createUI()
{
    return new FildesUI;
}

END_NAMESPACE_DISTRHO
