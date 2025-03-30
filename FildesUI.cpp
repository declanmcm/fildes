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
    }

    void setDenominator(const std::vector<double>& coeffs) {
        mDenominator = coeffs;
        mNeedsUpdate = true;
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
            
            // Rebuild the coefficients from roots
            std::vector<double>& coeffs = isPole ? mDenominator : mNumerator;
            
            // Convert roots back to coefficients by expanding (z - r_1)(z - r_2)...
            // For demonstration, we'll just update the relevant coefficient
            // In a real implementation, multiply out the polynomial terms
            
            // TODO: Implement proper roots-to-coefficients conversion
            // For now, just mark for update
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
                // Term for z^j * z
                newCoeffs[j + 1] += coeffs[j];
                
                // Term for z^j * (-r_i)
                newCoeffs[j] -= coeffs[j] * roots[i];
            }
            
            coeffs = newCoeffs;
        }
        
        // Convert complex coefficients to real (ignoring any imaginary part)
        std::vector<double> realCoeffs;
        for (const auto& c : coeffs) {
            realCoeffs.push_back(c.real());
        }
        
        return realCoeffs;
    }

    // Update coefficients from poles and zeros
    bool updateCoefficientsFromRoots() {
        // Convert zeros to numerator coefficients
        std::vector<double> newNum = rootsToCoeffs(mZeros);
        
        // Convert poles to denominator coefficients
        std::vector<double> newDen = rootsToCoeffs(mPoles);
        
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
        fCoeffT->setAbsolutePos(50, 50);
        fCoeffT->setSize(500, 30);
        fCoeffT->setCallback(this);

        fCoeffB = new EditableText(this, false);
        fCoeffB->setAbsolutePos(50, 100);
        fCoeffB->setSize(500, 30);
        fCoeffB->setCallback(this);

        // Set up the response and pole-zero areas
        fResponseArea = DGL::Rectangle<int>(50, 150, 700, 200);
        fPoleZeroArea = DGL::Rectangle<int>(500, 370, 250, 250);
        
        // Initialize dragging state
        fIsDragging = false;
        fDraggingPole = false;
        fDraggedIndex = -1;

        // we can use this if/when our resources are scalable, for now they are PNGs
        const double scaleFactor = getScaleFactor();
        setSize(800,400);
        if (scaleFactor != 1.0)
            setSize(400 * scaleFactor, 800 * scaleFactor);
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
        cairo_show_text(cr, "Transfer Function:"); // Static text

        drawFrequencyResponse(cr);
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
    
    void parseCoefficients(const char* str, std::vector<double>& coeffs) {
        coeffs.clear();
        
        // Create a copy of the input string to tokenize
        char inputCopy[1024]; // Ensure this is large enough
        strncpy(inputCopy, str, sizeof(inputCopy) - 1);
        inputCopy[sizeof(inputCopy) - 1] = '\0'; // Ensure null-termination
        
        // Tokenize the string using commas as delimiter
        char* token = strtok(inputCopy, ",");
        while (token != NULL) {
            // Convert token to double, treating empty tokens as 0
            double value = (token[0] == '\0') ? 0.0 : std::strtod(token, NULL);
            coeffs.push_back(value);
            
            // Move to the next token
            token = strtok(NULL, ",");
        }
    }

    void updateTransferFunction() {
        // Update the visualizer with current coefficients
        fVisualizer.setNumerator(fNumerator);
        fVisualizer.setDenominator(fDenominator);
        
        // Force a repaint to update the graph
        repaint();
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
            parseCoefficients(coeffStr, fNumerator);
            updateTransferFunction();
            setState("topTF", &newText[1]);
        }
        else if (newText[0] == 'B') {
            const char* coeffStr = &newText[1];
            parseCoefficients(coeffStr, fDenominator);
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
