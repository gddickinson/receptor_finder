#!/usr/bin/env python3
"""
Quick launcher for Receptor Finder GUI.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

if __name__ == "__main__":
    try:
        # Check tkinter first
        try:
            import tkinter
        except ImportError:
            print("\n❌ ERROR: tkinter is not available")
            print("\ntkinter is required for the GUI but is not installed.")
            print("\nHow to fix:")
            print("  • macOS: Install Python from python.org (includes tkinter)")
            print("  • Linux: sudo apt-get install python3-tk")
            print("  • Windows: Reinstall Python with tkinter option checked")
            print("\nOr run: python check_dependencies.py")
            sys.exit(1)
        
        # Import and run GUI
        from gui_app import main
        main()
        
    except ImportError as e:
        print("\n❌ ERROR: Missing dependencies")
        print("\nPlease install required packages:")
        print("  pip install -r requirements.txt")
        print("\nOr check what's missing:")
        print("  python check_dependencies.py")
        print(f"\nDetails: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
