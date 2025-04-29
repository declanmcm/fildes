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
 #include "Eigen/Eigen"

 
// RBJ Cookbook Filter Generator class (to replace existing FilterGenerator)
class RBJFilterGenerator {
    public:
        enum FilterType {
            LOW_PASS,
            HIGH_PASS,
            BAND_PASS,
            NOTCH,
            PEAK,
            LOW_SHELF,
            HIGH_SHELF,
            ALL_PASS
        };
        
        RBJFilterGenerator() 
            : mSampleRate(44100.0) {}
        
        void setSampleRate(double sampleRate) {
            mSampleRate = sampleRate;
        }
        
        // Generate filter coefficients using RBJ cookbook formulas
        std::pair<std::vector<double>, std::vector<double>> generateFilter(
            FilterType filterType,
            double frequency,
            double Q = 0.7071, // Default to Butterworth Q
            double gainDB = 0.0,
            double bandwidth = 1.0) // Bandwidth in octaves for shelf filters
        {
            // Normalize frequency (0 to π)
            double omega = 2.0 * M_PI * frequency / mSampleRate;
            double cosw = cos(omega);
            double sinw = sin(omega);
            double alpha;
            
            // Calculate alpha based on filter type
            if (filterType == PEAK || filterType == LOW_SHELF || filterType == HIGH_SHELF) {
                // For shelving and peaking filters, alpha can use bandwidth
                alpha = sinw * sinh(log(2.0) / 2.0 * bandwidth * omega / sinw);
            } else {
                // For other filters, use Q
                alpha = sinw / (2.0 * Q);
            }
            
            // Convert gain from dB to linear
            double A = pow(10.0, gainDB / 40.0); // Convert gainDB to linear amplitude
            
            // Initialize coefficients
            double b0 = 0.0, b1 = 0.0, b2 = 0.0;
            double a0 = 1.0, a1 = 0.0, a2 = 0.0;
            
            // Calculate coefficients based on filter type
            switch (filterType) {
                case LOW_PASS:
                    b0 = (1.0 - cosw) / 2.0;
                    b1 = 1.0 - cosw;
                    b2 = (1.0 - cosw) / 2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha;
                    break;
                    
                case HIGH_PASS:
                    b0 = (1.0 + cosw) / 2.0;
                    b1 = -(1.0 + cosw);
                    b2 = (1.0 + cosw) / 2.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha;
                    break;
                    
                case BAND_PASS: // Constant skirt gain (peak gain = Q)
                    b0 = alpha;
                    b1 = 0.0;
                    b2 = -alpha;
                    a0 = 1.0 + alpha;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha;
                    break;
                    
                case NOTCH:
                    b0 = 1.0;
                    b1 = -2.0 * cosw;
                    b2 = 1.0;
                    a0 = 1.0 + alpha;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha;
                    break;
                    
                case PEAK:
                    b0 = 1.0 + alpha * A;
                    b1 = -2.0 * cosw;
                    b2 = 1.0 - alpha * A;
                    a0 = 1.0 + alpha / A;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha / A;
                    break;
                    
                case LOW_SHELF:
                    {
                        double sqrtA = sqrt(A);
                        b0 = A * ((A + 1.0) - (A - 1.0) * cosw + 2.0 * sqrtA * alpha);
                        b1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * cosw);
                        b2 = A * ((A + 1.0) - (A - 1.0) * cosw - 2.0 * sqrtA * alpha);
                        a0 = (A + 1.0) + (A - 1.0) * cosw + 2.0 * sqrtA * alpha;
                        a1 = -2.0 * ((A - 1.0) + (A + 1.0) * cosw);
                        a2 = (A + 1.0) + (A - 1.0) * cosw - 2.0 * sqrtA * alpha;
                    }
                    break;
                    
                case HIGH_SHELF:
                    {
                        double sqrtA = sqrt(A);
                        b0 = A * ((A + 1.0) + (A - 1.0) * cosw + 2.0 * sqrtA * alpha);
                        b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cosw);
                        b2 = A * ((A + 1.0) + (A - 1.0) * cosw - 2.0 * sqrtA * alpha);
                        a0 = (A + 1.0) - (A - 1.0) * cosw + 2.0 * sqrtA * alpha;
                        a1 = 2.0 * ((A - 1.0) - (A + 1.0) * cosw);
                        a2 = (A + 1.0) - (A - 1.0) * cosw - 2.0 * sqrtA * alpha;
                    }
                    break;
                    
                case ALL_PASS:
                    b0 = 1.0 - alpha;
                    b1 = -2.0 * cosw;
                    b2 = 1.0 + alpha;
                    a0 = 1.0 + alpha;
                    a1 = -2.0 * cosw;
                    a2 = 1.0 - alpha;
                    break;
            }
            
            // Normalize coefficients by a0
            std::vector<double> numerator = {b0 / a0, b1 / a0, b2 / a0};
            std::vector<double> denominator = {1.0, a1 / a0, a2 / a0};
            
            // Debug output
            std::cout << "Generated RBJ filter coefficients:" << std::endl;
            std::cout << "Numerator: ";
            for (double c : numerator) std::cout << c << " ";
            std::cout << std::endl;
            std::cout << "Denominator: ";
            for (double c : denominator) std::cout << c << " ";
            std::cout << std::endl;
            
            return std::make_pair(numerator, denominator);
        }
        
    private:
        double mSampleRate;
    };
    
    // Add this class for a custom button widget (to be used for the Generate button)
class CustomButton : public CairoSubWidget,
                        public ButtonEventHandler,
                        public ButtonEventHandler::Callback
    {
    public:
        class Callback
        {
        public:
            virtual ~Callback() {}
            virtual void buttonClicked(CustomButton* button) = 0;
        };
    
        explicit CustomButton(SubWidget* const parent, const char* label)
            : CairoSubWidget(parent),
              ButtonEventHandler(this),
              fCallback(nullptr),
              fIsHovered(false),
              fIsPressed(false),
              fLabel(label)
        {
            ButtonEventHandler::setCallback(this);
        }

        explicit CustomButton(TopLevelWidget* const parent, const char* label)
            : CairoSubWidget(parent),
              ButtonEventHandler(this),
              fCallback(nullptr),
              fIsHovered(false),
              fIsPressed(false),
              fLabel(label)
        {
            ButtonEventHandler::setCallback(this);
        }
    
        void setCallback(Callback* const callback) noexcept
        {
            fCallback = callback;
        }

        void setId(uint32_t id) noexcept
        {
            fId = id;
        }

        uint32_t getId() const noexcept
        {
            return fId;
        }
    
    protected:
        void onCairoDisplay(const CairoGraphicsContext& context) override
        {
            cairo_t* const cr = context.handle;
    
            // Draw button background with different colors based on state
            if (fIsPressed) {
                cairo_set_source_rgb(cr, 0.4, 0.4, 0.8); // Darker blue when pressed
            } else if (fIsHovered) {
                cairo_set_source_rgb(cr, 0.5, 0.5, 0.9); // Lighter blue when hovered
            } else {
                cairo_set_source_rgb(cr, 0.6, 0.6, 1.0); // Default blue
            }
            
            // Draw rounded rectangle
            double radius = 10.0;
            double x = 0;
            double y = 0;
            double width = getWidth();
            double height = getHeight();
            
            // Rounded corners
            cairo_new_path(cr);
            cairo_arc(cr, x + width - radius, y + radius, radius, -M_PI/2, 0);
            cairo_arc(cr, x + width - radius, y + height - radius, radius, 0, M_PI/2);
            cairo_arc(cr, x + radius, y + height - radius, radius, M_PI/2, M_PI);
            cairo_arc(cr, x + radius, y + radius, radius, M_PI, 3*M_PI/2);
            cairo_close_path(cr);
            
            cairo_fill_preserve(cr);
            
            // Draw border
            cairo_set_source_rgb(cr, 0.2, 0.2, 0.7);
            cairo_set_line_width(cr, 2.0);
            cairo_stroke(cr);
    
            // Draw button text
            cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); // White text
            cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, 16.0);
            
            // Center text
            cairo_text_extents_t extents;
            cairo_text_extents(cr, fLabel.c_str(), &extents);
            
            double textX = (width - extents.width) / 2 - extents.x_bearing;
            double textY = (height - extents.height) / 2 - extents.y_bearing;
            
            cairo_move_to(cr, textX, textY);
            cairo_show_text(cr, fLabel.c_str());
        }
    
        bool onMouse(const MouseEvent& ev) override
        {
            return ButtonEventHandler::mouseEvent(ev);
        }
        
        // Fix: Use proper signature for buttonClicked
        void buttonClicked(SubWidget*, int button) override
        {
            if (button != 1) return;  // Only handle left-click
            
            if (fCallback != nullptr) {
                fCallback->buttonClicked(this);
            }
        }
        
        bool onMotion(const MotionEvent& ev) override
        {
            // Update hover state based on motion
            const bool wasHovered = fIsHovered;
            
            // Check if mouse is inside widget bounds
            fIsHovered = ev.pos.getX() >= 0 && ev.pos.getX() < getWidth() &&
                        ev.pos.getY() >= 0 && ev.pos.getY() < getHeight();
            
            // Repaint only if hover state changed
            if (fIsHovered != wasHovered) {
                repaint();
            }
            
            return CairoSubWidget::onMotion(ev);
        }
        
    private:
        Callback* fCallback;
        bool fIsHovered;
        bool fIsPressed;
        std::string fLabel;
        uint32_t fId;
    };

class Dropdown : public CairoSubWidget,
                public ButtonEventHandler,
                public ButtonEventHandler::Callback
    {
    public:
        class Callback
        {
        public:
            virtual ~Callback() {}
            virtual void selectionChanged(Dropdown* dropdown, int selectedIndex) = 0;
        };

        explicit Dropdown(SubWidget* const parent, const std::vector<std::string>& items)
            : CairoSubWidget(parent),
                ButtonEventHandler(this),
                fCallback(nullptr),
                fItems(items),
                fSelectedIndex(0),
                fIsOpen(false),
                fIsHovered(false),
                fClosedHeight(30)
        {
            ButtonEventHandler::setCallback(this);
            setSize(getWidth(), fClosedHeight);
        }

        explicit Dropdown(TopLevelWidget* const parent, const std::vector<std::string>& items)
            : CairoSubWidget(parent),
                ButtonEventHandler(this),
                fCallback(nullptr),
                fItems(items),
                fSelectedIndex(0),
                fIsOpen(false),
                fIsHovered(false),
                fClosedHeight(30)
        {
            ButtonEventHandler::setCallback(this);
            setSize(getWidth(), fClosedHeight);
        }

        int getSelectedIndex() const
        {
            return fSelectedIndex;
        }

        void setSelectedIndex(int index)
        {
            if (index >= 0 && index < static_cast<int>(fItems.size())) {
                fSelectedIndex = index;
                repaint();
            }
        }

        const std::string& getSelectedItem() const
        {
            return fItems[fSelectedIndex];
        }

        void setCallback(Callback* const callback) noexcept
        {
            fCallback = callback;
        }

        void setId(uint32_t id) noexcept
        {
            fId = id;
        }

        uint32_t getId() const noexcept
        {
            return fId;
        }

        // Override setSize to store width for resizing during open/close
        void setSize(const uint width, const uint height)
        {
            fWidth = width;
            CairoSubWidget::setSize(width, height);
        }

        // Override setAbsolutePos to store position for dropdown positioning
        void setAbsolutePos(const int x, const int y)
        {
            fPosX = x;
            fPosY = y;
            CairoSubWidget::setAbsolutePos(x, y);
        }

    protected:
        void onCairoDisplay(const CairoGraphicsContext& context) override
        {
            cairo_t* const cr = context.handle;
            
            // Draw dropdown box
            if (fIsHovered) {
                cairo_set_source_rgb(cr, 0.95, 0.95, 0.95); // Lighter when hovered
            } else {
                cairo_set_source_rgb(cr, 0.9, 0.9, 0.9); // Default gray
            }
            
            // Main box
            cairo_rectangle(cr, 0, 0, getWidth(), 30);
            cairo_fill_preserve(cr);
            
            // Border
            cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
            cairo_set_line_width(cr, 1.0);
            cairo_stroke(cr);
            
            // Draw selected item text
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
            cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
            cairo_set_font_size(cr, 14.0);
            cairo_move_to(cr, 10, 20);
            cairo_show_text(cr, fItems[fSelectedIndex].c_str());
            
            // Draw dropdown arrow
            cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
            cairo_move_to(cr, getWidth() - 20, 10);
            cairo_line_to(cr, getWidth() - 10, 10);
            cairo_line_to(cr, getWidth() - 15, 20);
            cairo_close_path(cr);
            cairo_fill(cr);
            
            // Draw dropdown menu if open
            if (fIsOpen) {
                const int itemHeight = 25;
                const int menuHeight = fItems.size() * itemHeight;
                
                // Menu background
                cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
                cairo_rectangle(cr, 0, 30, getWidth(), menuHeight);
                cairo_fill_preserve(cr);
                
                // Menu border
                cairo_set_source_rgb(cr, 0.4, 0.4, 0.4);
                cairo_stroke(cr);
                
                // Draw items
                for (size_t i = 0; i < fItems.size(); i++) {
                    const int y = 30 + i * itemHeight;
                    
                    // Highlight selected item
                    if (static_cast<int>(i) == fSelectedIndex) {
                        cairo_set_source_rgb(cr, 0.8, 0.8, 0.9);
                        cairo_rectangle(cr, 0, y, getWidth(), itemHeight);
                        cairo_fill(cr);
                    }
                    
                    // Item text
                    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
                    cairo_move_to(cr, 10, y + 18);
                    cairo_show_text(cr, fItems[i].c_str());
                    
                    // Separator line
                    if (i < fItems.size() - 1) {
                        cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
                        cairo_move_to(cr, 0, y + itemHeight);
                        cairo_line_to(cr, getWidth(), y + itemHeight);
                        cairo_stroke(cr);
                    }
                }
            }
        }

        bool onMouse(const MouseEvent& ev) override
        {
            const int y = ev.pos.getY();
            
            if (fIsOpen && y > 30) {
                // Click in dropdown menu
                if (ev.press && ev.button == 1) {
                    const int itemHeight = 25;
                    const int itemIndex = (y - 30) / itemHeight;
                    
                    if (itemIndex >= 0 && itemIndex < static_cast<int>(fItems.size())) {
                        // Select item

                        toggleDropdown(false);
                        fSelectedIndex = itemIndex;
                        repaint();
                        
                        if (fCallback) {
                            fCallback->selectionChanged(this, fSelectedIndex);
                        }
                    }
                    return true;
                }
            } else {
                return ButtonEventHandler::mouseEvent(ev);
            }
            
            return false;
        }
        
        void buttonClicked(SubWidget* widget, int button) override
        {
            if (button == 1) {
                // Toggle dropdown
                toggleDropdown(!fIsOpen);
            }
        }
        
        bool onMotion(const MotionEvent& ev) override
        {
            const int y = ev.pos.getY();
            
            // Update hover state
            const bool wasHovered = fIsHovered;
            fIsHovered = ev.pos.getX() >= 0 && ev.pos.getX() < getWidth() &&
                        ev.pos.getY() >= 0 && ev.pos.getY() < getHeight();
                        
            if (wasHovered != fIsHovered) {
                repaint();
            }
            
            if (fIsOpen && y > 30) {
                // Hover in dropdown menu
                const int itemHeight = 25;
                const int hoverIndex = (y - 30) / itemHeight;
                
                if (hoverIndex >= 0 && hoverIndex < static_cast<int>(fItems.size())) {
                    // Highlight hovered item (simplified for now)
                    repaint();
                }
            }
            
            return CairoSubWidget::onMotion(ev);
        }

    private:

        void toggleDropdown(bool open) {            
            fIsOpen = open;
            
            const int itemHeight = 25;
            const int menuHeight = fItems.size() * itemHeight;
            
            if (open) {
                toFront();
                // Expand height to show menu
                setSize(fWidth, fClosedHeight + menuHeight);
            } else {
                // Shrink back to closed height
                setSize(fWidth, fClosedHeight);
            }
            
            repaint();
        }

        Callback* fCallback;
        std::vector<std::string> fItems;
        int fSelectedIndex;
        bool fIsOpen;
        bool fIsHovered;
        uint32_t fId;
        uint fWidth;
        int fPosX;
        int fPosY;
        const int fClosedHeight;
    };

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
                 roots.push_back(std::complex<double>(-c[1] / c[0], 0.0));
             }
             return roots;
         }
         
         // Normalize the coefficients (a_n = 1)
         double lead = c.back();
         for (double& coef : c) {
             coef /= lead;
         }
         
         const int n = c.size() - 1;
         
         // Create companion matrix
         Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);  
         
         // Fill the first row with negated coefficients
         for (int i = 0; i < n; i++) {
             A(n - 1 - i, n - 1) = -c[i + 1] / c[0];
         }
         
         // Fill the subdiagonal with ones
         for (int i = 1; i < n; i++) {
             A(i, i - 1) = 1.0;
         }
 
         std::cout << "matrix:" << A << std::endl;
         
         // Compute eigenvalues (roots)
         Eigen::EigenSolver<Eigen::MatrixXd> es(A);
         
         // Extract the eigenvalues
         for (int i = 0; i < n; ++i) {
             std::complex<double> root(es.eigenvalues()[i].real(), es.eigenvalues()[i].imag());
             roots.push_back(root);
         }
         
         return roots;
     }
 
     void findZeros() {
         // Find zeros of the transfer function (roots of numerator)
         mZeros = findRoots(mNumerator);
         updateConjugateIndices(mZeros, mZeroConjugateIndices);
     }
 
     void findPoles() {
         // Find poles of the transfer function (roots of denominator)
         mPoles = findRoots(mDenominator);
         updateConjugateIndices(mPoles, mPoleConjugateIndices);
     }
 
     void updateConjugateIndices(const std::vector<std::complex<double>>& roots, 
         std::vector<int>& conjugateIndices) {
         // Reset conjugate indices array
         const size_t numRoots = roots.size();
         conjugateIndices.assign(numRoots, -1);
 
         // Find conjugate pairs
         for (size_t i = 0; i < numRoots; ++i) {
             // Skip if this root already has a conjugate assigned
             if (conjugateIndices[i] != -1) continue;
 
             // Check if this root has a significant imaginary part
             if (std::abs(roots[i].imag()) > 1e-10) {
                 // Look for its conjugate
                 for (size_t j = i + 1; j < numRoots; ++j) {
                     double realDiff = std::abs(roots[j].real() - roots[i].real());
                     double imagSum = std::abs(roots[j].imag() + roots[i].imag());
 
                     if (realDiff < 1e-10 && imagSum < 1e-10) {
                         // Found conjugate pair, store bidirectional indices
                         conjugateIndices[i] = j;
                         conjugateIndices[j] = i;
                         break;
                     }
                 }
             }
         }
 
         // Debug output
         std::cout << "Conjugate Indices: ";
         for (int idx : conjugateIndices) {
         std::cout << idx << " ";
         }
         std::cout << std::endl;
     }
 
     // Update a pole or zero based on dragging
     void updatePoleZero(bool isPole, int index, const std::complex<double>& newValue) {
         if (index < 0) return;
         
         // Get references to the appropriate arrays
         std::vector<std::complex<double>>& roots = isPole ? mPoles : mZeros;
         std::vector<int>& conjugateIndices = isPole ? mPoleConjugateIndices : mZeroConjugateIndices;
         
         if (index >= static_cast<int>(roots.size())) return;
         
         // Update the root
         roots[index] = newValue;
         
         // Check if this root has a conjugate pair that needs updating
         int conjugateIndex = conjugateIndices[index];
         if (conjugateIndex >= 0 && conjugateIndex < static_cast<int>(roots.size())) {
             // Update the conjugate with the appropriate value
             roots[conjugateIndex] = std::complex<double>(newValue.real(), -newValue.imag());
             
             // Debug output
             std::cout << "Updated conjugate at index " << conjugateIndex << std::endl;
         } else if (std::abs(newValue.imag()) > 1e-10) {
             // This root has a significant imaginary part, but no conjugate is assigned
             // Try to find a real root to convert to a conjugate, or the closest root
             int newConjugateIndex = -1;
             double minDistance = 1e10;
             
             for (size_t i = 0; i < roots.size(); ++i) {
                 if (i != index && conjugateIndices[i] == -1) {
                     // Prioritize real roots (imaginary part close to zero)
                     if (std::abs(roots[i].imag()) < 1e-10) {
                         newConjugateIndex = i;
                         break;
                     }
                     
                     // Otherwise find the closest root
                     double distance = std::abs(roots[i] - std::complex<double>(newValue.real(), -newValue.imag()));
                     if (distance < minDistance) {
                         minDistance = distance;
                         newConjugateIndex = i;
                     }
                 }
             }
             
             // If we found a suitable root, convert it to a conjugate
             if (newConjugateIndex >= 0) {
                 roots[newConjugateIndex] = std::complex<double>(newValue.real(), -newValue.imag());
                 
                 // Update conjugate indices in both directions
                 conjugateIndices[index] = newConjugateIndex;
                 conjugateIndices[newConjugateIndex] = index;
                 
                 std::cout << "Created new conjugate pair: " << index << " <-> " << newConjugateIndex << std::endl;
             } else {
                 // Couldn't find a suitable conjugate, force this to be real
                 roots[index] = std::complex<double>(newValue.real(), 0.0);
                 conjugateIndices[index] = -1;
                 
                 std::cout << "Forced root to be real - no conjugate available" << std::endl;
             }
         }
 
         std::cout << "roots: ";
         for (const auto& coeff : roots) {
             std::cout << coeff << " ";
         }
         std::cout << std::endl;
         
         // Immediately update the coefficients
         if (isPole) {
             std::vector<double> newDenominator = rootsToCoeffs(mPoles);
             std::cout << "new den: ";
             for (const auto& coeff : newDenominator) {
                 std::cout << coeff << " ";
             }
             std::cout << std::endl;
             
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
             std::cout << "new num: ";
             for (const auto& coeff : mNumerator) {
                 std::cout << coeff << " ";
             }
             std::cout << std::endl;
         }
         
         // Mark for response update
         mNeedsUpdate = true;
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
 
 
         std::reverse(realCoeffs.begin(), realCoeffs.end());
         
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
 
     int getConjugateIndex(bool isPole, int index) const {
         if (index < 0) return -1;
         
         const std::vector<int>& indices = isPole ? mPoleConjugateIndices : mZeroConjugateIndices;
         
         if (index < static_cast<int>(indices.size())) {
             return indices[index];
         }
         
         return -1;
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
     std::vector<int> mZeroConjugateIndices;
     std::vector<int> mPoleConjugateIndices;
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
                        public EditableText::Callback,
                        public CustomButton::Callback,
                        public Dropdown::Callback
 {
     ScopedPointer<CairoImageKnob> fKnob;
     ScopedPointer<CairoImageSwitch> fButton;
     ScopedPointer<DemoWidgetClickable> fWidgetClickable;
     ScopedPointer<EditableText> fCoeffT, fCoeffB;
     ScopedPointer<DGL::SubWidget> fContainer;

     ScopedPointer<RBJFilterGenerator> fRBJFilterGenerator;
     ScopedPointer<Dropdown> fFilterTypeDropdown;
     ScopedPointer<EditableText> fFrequencyInput;
     ScopedPointer<EditableText> fQInput;
     ScopedPointer<EditableText> fGainDBInput;
     ScopedPointer<EditableText> fBandwidthInput;
     ScopedPointer<CustomButton> fGenerateButton;


    ScopedPointer<EditableText> *fTextInputs[7];
 
     // Add the visualizer
     TransferFunctionVisualizer fVisualizer;

     ScopedPointer<EditableText> fGainInput;
     
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
     bool fGeneratorTracking;
     double fDragStartX, fDragStartY;
 
     int fDraggedIndex;
 
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
         fGeneratorTracking = false;
         
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
         fCoeffT = new EditableText(this, "T");
         fCoeffT->setAbsolutePos(150, 50);
         fCoeffT->setSize(500, 30);
         fCoeffT->setCallback(this);
 
         fCoeffB = new EditableText(this, "B");
         fCoeffB->setAbsolutePos(150, 100);
         fCoeffB->setSize(500, 30);
         fCoeffB->setCallback(this);

         // Create gain input field
        fGainInput = new EditableText(this, "G");
        fGainInput->setAbsolutePos(670, 100);
        fGainInput->setSize(70, 30);
        fGainInput->setCallback(this);
        fGainInput->setText("");
         
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

        // Initialize RBJ filter generator
        fRBJFilterGenerator = new RBJFilterGenerator();
        fRBJFilterGenerator->setSampleRate(getSampleRate());

        // Create filter type dropdown for RBJ cookbook filters
        std::vector<std::string> filterTypes = {
            "Low Pass", "High Pass", "Band Pass", "Notch", 
            "Peak", "Low Shelf", "High Shelf", "All Pass"
        };
        fFilterTypeDropdown = new Dropdown(this, filterTypes);
        fFilterTypeDropdown->setAbsolutePos(150, 420);
        fFilterTypeDropdown->setSize(200, 30);
        fFilterTypeDropdown->setCallback(this);

        // Create frequency input
        fFrequencyInput = new EditableText(this, " ");
        fFrequencyInput->setAbsolutePos(150, 470);
        fFrequencyInput->setSize(150, 30);
        fFrequencyInput->setCallback(this);
        fFrequencyInput->setText("1000");  // Default value

        // Create Q input
        fQInput = new EditableText(this, " ");
        fQInput->setAbsolutePos(150, 520);
        fQInput->setSize(150, 30);
        fQInput->setCallback(this);
        fQInput->setText("0.7071");  // Default Butterworth Q

        // Create gain dB input (for peak and shelf filters)
        fGainDBInput = new EditableText(this, " ");
        fGainDBInput->setAbsolutePos(150, 570);
        fGainDBInput->setSize(150, 30);
        fGainDBInput->setCallback(this);
        fGainDBInput->setText("0.0");  // Default value
        fGainDBInput->setVisible(false);

        // Create bandwidth input (alternative to Q for some filters)
        fBandwidthInput = new EditableText(this, " ");
        fBandwidthInput->setAbsolutePos(320, 520);
        fBandwidthInput->setSize(150, 30);
        fBandwidthInput->setCallback(this);
        fBandwidthInput->setText("1.0");  // Default value (in octaves)
        fBandwidthInput->setVisible(false);

        // Create generate button
        fGenerateButton = new CustomButton(this, "Toggle Tracking");
        fGenerateButton->setAbsolutePos(150, 620);
        fGenerateButton->setSize(150, 40);
        fGenerateButton->setCallback(this);

        // Add these to the fTextInputs array
        fTextInputs[3] = &fFrequencyInput;
        fTextInputs[4] = &fQInput;
        fTextInputs[5] = &fGainDBInput;
        fTextInputs[6] = &fBandwidthInput;
        fTextInputs[0] = &fCoeffT;
        fTextInputs[1] = &fCoeffB;
        fTextInputs[2] = &fGainInput;

        // Initialize visibility based on current filter type
        updateFilterParametersVisibility();
     }
 
 protected:

    void drawFilterParameterLabels(cairo_t* cr)
    {
        // Draw filter generator section title
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        if (fGeneratorTracking) cairo_set_source_rgb(cr, 0.0, 0.0, 0.5);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 20);
        cairo_move_to(cr, 10, 400);
        cairo_show_text(cr, "Filter Generator");
        
        // Draw filter parameter labels
        cairo_set_font_size(cr, 16);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        
        cairo_move_to(cr, 10, 440);
        cairo_show_text(cr, "Filter Type:");
        
        cairo_move_to(cr, 10, 490);
        cairo_show_text(cr, "Frequency (Hz):");
        
        cairo_move_to(cr, 10, 540);
        cairo_show_text(cr, "Q:");
        
        // Show bandwidth label if applicable
        if (fBandwidthInput->isVisible()) {
            cairo_move_to(cr, 320, 490);
            cairo_show_text(cr, "Bandwidth (oct):");
        }
        
        // Show gain label if applicable
        if (fGainDBInput->isVisible()) {
            cairo_move_to(cr, 10, 540);
            cairo_show_text(cr, "Gain (dB):");
        }
    }

     void onCairoDisplay(const CairoGraphicsContext& context)
     {
         cairo_t* const cr = context.handle;
         cairo_set_source_rgb(cr, 0.541, 0.576, 1);
         cairo_paint(cr);
 
         cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); // Set text color
         cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
         cairo_set_font_size(cr, 20);
 
         cairo_move_to(cr, 10, 30);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 20);
        cairo_show_text(cr, "Transfer Function");
         
         cairo_move_to(cr, 10, 75);
         cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
         cairo_show_text(cr, "Numerator:");
         
         cairo_move_to(cr, 10, 125);
         cairo_show_text(cr, "Denominator:");
 
         drawFrequencyResponse(cr);
 
         drawPoleZeroPlot(cr);
        drawFilterParameterLabels(cr);
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
             // For denominator, ALWAYS add z^0 term with coefficient 1.0 at the beginning
             // This is the implicit term that's not displayed in the text field
             coeffs.push_back(1.0);
             
             // Then add the parsed coefficients (which are z^-1 and higher terms)
             coeffs.insert(coeffs.end(), tempCoeffs.begin(), tempCoeffs.end());
         } else {
             // For numerator, use the coefficients as is
             coeffs = tempCoeffs;
         }
         
         // Debug: print the parsed coefficients
         std::cout << "Parsed coefficients (isDenominator=" << isDenominator << "): ";
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
                 index = static_cast<int>(i);
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
                 index = static_cast<int>(i);
             }
         }
         
         return closest < minDistance;
     }
     
     // Convert screen coordinates to complex plane coordinates
     std::complex<double> screenToComplex(int x, int y, bool maintainConjugatePairs = true) {
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
         
         // Debug output
         std::cout << "Drag to: (" << newValue.real() << ", " << newValue.imag() << ")" << std::endl;
         
         // Update pole or zero using visualizer's method which handles conjugate pairs
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
 
     int findConjugateIndex(bool isPole, int index) {
         const std::vector<std::complex<double>>& roots = 
             isPole ? fVisualizer.getPoles() : fVisualizer.getZeros();
         
         if (index >= 0 && index < static_cast<int>(roots.size()) && 
             std::abs(roots[index].imag()) > 1e-6) {
             for (int i = 0; i < static_cast<int>(roots.size()); ++i) {
                 if (i != index) {
                     double realDiff = std::abs(roots[i].real() - roots[index].real());
                     double imagSum = std::abs(roots[i].imag() + roots[index].imag());
                     
                     if (realDiff < 1e-6 && imagSum < 1e-6) {
                         return i;
                     }
                 }
             }
         }
         
         return -1;
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
                         
                         std::cout << "Drag START: " << (isPole ? "POLE" : "ZERO") 
                                   << " index " << index << std::endl;
                         
                         return true;
                     }
                 } else {
                     // Mouse up
                     if (fIsDragging) {
                         std::cout << "Drag END" << std::endl;
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
        for (int i = 0; i < 7; i++) {
            if (*(fTextInputs[i]) != nullptr) {
                (*(fTextInputs[i]))->setFocused(*(fTextInputs[i]) == current);
            }
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

     // Implement CustomButton::Callback
     void buttonClicked(CustomButton* button)
    {
        if (button == fGenerateButton.get()) {
            generateFilter();
            fGeneratorTracking = true;
        }
    }
    
    // Dropdown::Callback
    void selectionChanged(Dropdown* dropdown, int selectedIndex) override
    {
        if (dropdown == fFilterTypeDropdown.get()) {
            // Update parameter visibility based on filter type
            updateFilterParametersVisibility();
            generateFilter();
            repaint();
        }
    }

     void setGain(std::string gainStr) override
    {
        // Get the gain value from the input field
        double gain = 1.0;
        
        try {
            gain = std::stod(gainStr);
            
            // Validate gain (prevent zero/negative values)
            if (gain <= 0.0) {
                std::cout << "Invalid gain value: " << gain << std::endl;
                fGainInput->setText(""); // Reset to default gain
                return;
            }
        } catch (const std::exception& e) {
            std::cout << "Error parsing gain: " << e.what() << std::endl;
            fGainInput->setText(""); // Reset to default gain
            return;
        }
        
        // Apply the gain by multiplying all numerator coefficients
        for (double& coeff : fNumerator) {
            coeff *= gain;
        }
        
        // Update the visualizer with modified coefficients
        fVisualizer.setNumerator(fNumerator);
        
        // Format coefficients for display and update the text field
        std::string numStr = formatCoefficients(fNumerator, false);
        fCoeffT->setText(numStr);
        
        // Update the plugin state
        setState("topTF", numStr.c_str());
        
        // Reset gain to 1.0 after applying
        fGainInput->setText("");
        
        // Force repaint to update the display
        repaint();
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
 private:

  // Update visibility of parameters based on filter type
    void updateFilterParametersVisibility()
    {
        // All filter types use frequency and Q
        fFrequencyInput->setVisible(true);
        fQInput->setVisible(true);
        
        // Only certain filter types use gain
        RBJFilterGenerator::FilterType filterType = 
            static_cast<RBJFilterGenerator::FilterType>(fFilterTypeDropdown->getSelectedIndex());
        
        switch (filterType) {
            case RBJFilterGenerator::PEAK:
            case RBJFilterGenerator::LOW_SHELF:
            case RBJFilterGenerator::HIGH_SHELF:
                fGainDBInput->setVisible(true);
                fBandwidthInput->setVisible(true);
                break;
                
            default:
                fGainDBInput->setVisible(false);
                fBandwidthInput->setVisible(false);
                break;
        }
    }

    void generateIfTracking() override {
        if (fGeneratorTracking)
            generateFilter();
    }
    
    void generateFilter()
    {
        // Get selected filter type
        RBJFilterGenerator::FilterType filterType = 
            static_cast<RBJFilterGenerator::FilterType>(fFilterTypeDropdown->getSelectedIndex());
            
        // Parse parameters from input fields
        double frequency = 1000.0;  // Default
        double Q = 0.7071;  // Default
        double gainDB = 0.0;  // Default
        double bandwidth = 1.0;  // Default
        
        try {
            // Parse frequency
            frequency = std::stod(fFrequencyInput->getText());
            frequency = std::max(20.0, std::min(20000.0, frequency));  // Clamp to audible range
            
            // Parse Q
            Q = std::stod(fQInput->getText());
            Q = std::max(0.1, std::min(30.0, Q));  // Reasonable Q range
            
            // Parse gain for peak/shelf filters
            if (fGainDBInput->isVisible()) {
                gainDB = std::stod(fGainDBInput->getText());
                gainDB = std::max(-30.0, std::min(30.0, gainDB));  // Reasonable gain range
            }
            
            // Parse bandwidth for peak/shelf filters
            if (fBandwidthInput->isVisible()) {
                bandwidth = std::stod(fBandwidthInput->getText());
                bandwidth = std::max(0.1, std::min(4.0, bandwidth));  // Reasonable bandwidth range
            }
        } catch (const std::exception& e) {
            std::cout << "Error parsing filter parameters: " << e.what() << std::endl;
            return;
        }
        
        // Generate filter coefficients
        std::pair<std::vector<double>, std::vector<double>> coeffs = 
            fRBJFilterGenerator->generateFilter(filterType, frequency, Q, gainDB, bandwidth);
        
        // Update the filter coefficients
        fNumerator = coeffs.first;
        fDenominator = coeffs.second;
        
        // Format coefficients for display
        std::string numStr = formatCoefficients(fNumerator, false);
        std::string denStr = formatCoefficients(fDenominator, true);
        
        // Update coefficient text fields
        fCoeffT->setText(numStr);
        fCoeffB->setText(denStr);
        
        // Update visualizer and GUI
        updateVisualizationFromCoefficients();
        
        // Update plugin state
        setState("topTF", numStr.c_str());
        setState("bottomTF", denStr.c_str());
        
        std::cout << "RBJ filter generated successfully!" << std::endl;
    }
 };
 
 UI* createUI()
 {
     return new FildesUI;
 }
 
 END_NAMESPACE_DISTRHO