#!/bin/bash
# install_dependencies.sh - Install all required Python dependencies for plotting

# Function to display colorful messages
function echo_success() { echo -e "\033[0;32m✓ $1\033[0m"; }
function echo_info() { echo -e "\033[0;34mℹ $1\033[0m"; }
function echo_warning() { echo -e "\033[0;33m⚠ $1\033[0m"; }
function echo_error() { echo -e "\033[0;31m❌ $1\033[0m"; }

echo_info "Checking and installing required Python dependencies for plotting..."

# Check if venv is active
if [[ -z "$VIRTUAL_ENV" ]]; then
    echo_warning "Virtual environment not active. Activating venv..."
    
    # Check if venv exists
    if [ -d "venv" ]; then
        source venv/bin/activate
        echo_success "Virtual environment activated"
    else
        echo_info "Creating new virtual environment..."
        python3 -m venv venv
        source venv/bin/activate
        echo_success "Virtual environment created and activated"
    fi
else
    echo_success "Virtual environment is already active: $VIRTUAL_ENV"
fi

# Update pip
echo_info "Updating pip..."
pip install --upgrade pip

# Install dependencies with progress output
function install_package() {
    echo_info "Installing $1..."
    pip install $1
    if [ $? -eq 0 ]; then
        echo_success "$1 installed successfully"
        return 0
    else
        echo_error "Failed to install $1"
        return 1
    fi
}

# Required packages for plotting
REQUIRED_PACKAGES=(
    "numpy"
    "matplotlib"
    "scipy"
    "pandas"
    "cycler"
    "kiwisolver"
    "Pillow"
    "pyparsing"
    "python-dateutil"
)

# Install packages
echo_info "Installing required Python packages..."
for package in "${REQUIRED_PACKAGES[@]}"; do
    install_package "$package"
done

# Check for FFmpeg for animations
echo_info "Checking for FFmpeg (required for animations)..."
if command -v ffmpeg &> /dev/null; then
    echo_success "FFmpeg found: $(which ffmpeg)"
else
    echo_warning "FFmpeg not found, which is required for animations."
    echo_info "You should install FFmpeg to enable animation generation:"
    
    # Different instructions based on OS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo_info "For macOS: brew install ffmpeg"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        if command -v apt-get &> /dev/null; then
            echo_info "For Ubuntu/Debian: sudo apt-get install ffmpeg"
        elif command -v yum &> /dev/null; then
            echo_info "For RHEL/CentOS: sudo yum install ffmpeg"
        fi
    fi
fi

# Test imports to verify installation
echo_info "Verifying Python package installation..."

python3 -c "
import sys
import numpy
import matplotlib
import scipy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.cm import get_cmap
import argparse
from pathlib import Path
import re

print('Python version:', sys.version)
print('NumPy version:', numpy.__version__)
print('Matplotlib version:', matplotlib.__version__)
print('SciPy version:', scipy.__version__)
print('All required modules imported successfully!')
"

if [ $? -eq 0 ]; then
    echo_success "All required packages are installed correctly."
else
    echo_error "There was a problem with the package installation or imports."
    echo_info "Please check the error messages above."
fi

echo_info "Testing plotting script imports..."
if [ -f "plot_pulsatile_flow.py" ]; then
    python3 -c "
import importlib.util
import sys

try:
    spec = importlib.util.spec_from_file_location('plot_module', 'plot_pulsatile_flow.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules['plot_module'] = module
    spec.loader.exec_module(module)
    print('Successfully imported plot_pulsatile_flow.py')
except Exception as e:
    print(f'Error importing plot_pulsatile_flow.py: {e}')
    exit(1)
"
    if [ $? -eq 0 ]; then
        echo_success "plot_pulsatile_flow.py imports correctly."
    else
        echo_error "There was a problem importing plot_pulsatile_flow.py."
    fi
else
    echo_warning "plot_pulsatile_flow.py not found in current directory."
fi

echo_info "All dependencies are now installed and verified!"