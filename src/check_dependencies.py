#!/usr/bin/env python3
"""
Check if all dependencies are installed and working.
"""
import sys

print("Checking dependencies for Receptor Finder GUI...\n")

# Check Python version
print(f"Python version: {sys.version}")
if sys.version_info < (3, 7):
    print("✗ Python 3.7+ required")
    sys.exit(1)
else:
    print("✓ Python version OK")

# Check required packages
required = {
    'numpy': 'NumPy',
    'scipy': 'SciPy', 
    'pandas': 'Pandas',
    'matplotlib': 'Matplotlib',
    'seaborn': 'Seaborn',
    'tkinter': 'Tkinter (GUI)',
}

missing = []
for module, name in required.items():
    try:
        __import__(module)
        print(f"✓ {name}")
    except ImportError:
        print(f"✗ {name} - NOT FOUND")
        missing.append(module)

if missing:
    print(f"\n❌ Missing packages: {', '.join(missing)}")
    print("\nTo install:")
    if 'tkinter' in missing:
        print("  For tkinter:")
        print("    macOS: Should be included with Python from python.org")
        print("    Linux: sudo apt-get install python3-tk")
        print("    Windows: Should be included with Python installer")
    
    other_missing = [m for m in missing if m != 'tkinter']
    if other_missing:
        print(f"\n  pip install {' '.join(other_missing)}")
    print("\nOr install all at once:")
    print("  pip install -r requirements.txt")
else:
    print("\n✅ All dependencies installed!")
    print("\nYou can now run:")
    print("  python launch_gui.py")
