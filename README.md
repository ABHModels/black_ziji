# BlackZiji

BlackZiji is a general relativistic ray-tracing and X-ray reflection modeling toolkit. It enables fast simulation of photon trajectories around black holes and supports convolution with reflection models like XILLVER. The codebase includes C++ high-performance backends and Python bindings for easy use.

---

## 📦 Features

- Supports both standard Kerr and modified spacetime metrics
- Lamp-post corona geometry with relativistic ray tracing
- Full relativistic X-ray reflection spectrum calculation
- Iron line



---


## 🚀 Installation

### 1. Clone the repository

```bash
git clone https://github.com/ABHModels/black_ziji.git
cd black_ziji
```


### 2. Install system-level dependencies

These packages are required to compile the C++ backend, enable OpenMP support, and build the Python C extensions using Cython.

#### For **Ubuntu/Debian**:

```bash
sudo apt update
sudo apt install build-essential cmake libgsl-dev libomp-dev python3-dev python3-pip
```

#### For **CentOS/RHEL/Fedora**:

```bash
sudo dnf install gcc gcc-c++ cmake gsl-devel libomp-devel python3-devel python3-pip
```

#### For **macOS (Homebrew)**:

```bash
brew install gsl libomp gcc cmake
```

> Make sure Python 3 and `pip` are installed and accessible. Use a virtual environment (`venv`) to avoid system-wide changes.


### 3. Set up Python environment (optional but recommended)

```bash
python3 -m venv venv
source venv/bin/activate
```

### 4. Install required Python packages

```bash
pip install -r requirements.txt
```

> **Tested with Python 3.8+**

#### `requirements.txt` contents:

```
astropy>=5.2.2
Cython>=3.0.9
matplotlib>=3.7.5
numpy>=1.24.4
releash>=0.5.1
setuptools>=56.0.0
```

---

## 🛠️ Building the Project

Run the following from the root of the project:

```bash
cmake .
make all
make line_bin
make reb_bin
```



## 🧪 Running the Python Interface

After building, you can run the main script:

```bash
python3 main.py
```

If you encounter `ModuleNotFoundError`, try:

```bash
PYTHONPATH=bin python3 main.py
```

---

## 📁 Project Structure

```
black_ziji/
├── bin/              # Compiled .so libraries
├── conv_core/        # Spectral convolution, Cython Python bindings
├── data/             # Input/output datasets
├── external/         # xtensor, xtl dependencies
├── main.py           # Python main entry point
├── zijiray/          # Core C++ ray tracing engine
├── xillver/          # XILLVER integration
├── CMakeLists.txt    # CMake configuration
├── Makefile          # Optional Makefile-based build
├── requirements.txt  # Python dependencies
└── README.md         # You're here!
```

---

## ⚠️ Troubleshooting

**Missing `Python.h` during Cython build?**  
→ Install dev headers:  
- Ubuntu/Debian: `sudo apt install python3-dev`  
- CentOS/Fedora: `sudo dnf install python3-devel`  

**Missing `libomp`?**  
→ Ubuntu: `sudo apt install libomp-dev`  
→ Fedora: `sudo dnf install libomp-devel`  

**Can't import `ray_line` in Python?**  
→ Ensure `bin/conv/ray_line.so` exists and use `PYTHONPATH=bin`



---

## 📬 Contact

