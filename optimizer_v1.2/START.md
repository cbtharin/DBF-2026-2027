# How to start it

## 1. Check the folder looks like this

```
optimizerv_1.2/          <-- you will run the command from HERE
└── dbf2027/            <-- this folder must be named exactly "dbf2027"
    ├── __init__.py
    ├── gui.py
    └── ... 18 more .py files
```

The outer folder can be called anything. The inner one cannot.

## 2. Run this

```bash
python -m dbf2027.gui
```

On macOS or Linux, use `python3` instead of `python`.

That's it. The window that opens is the whole program.

## If it doesn't start

| what you see | what's wrong | fix |
|---|---|---|
| `No module named 'dbf2027'` | you're in the wrong folder, or the folder isn't named `dbf2027` | `cd` one level *above* `dbf2027`, and check the spelling |
| `attempted relative import with no known parent package` | you ran `python gui.py` | use the command in step 2 — `gui.py` can't run on its own |
| `No module named 'tkinter'` | Python was installed without tkinter | Linux: `sudo apt install python3-tk`. Windows/macOS: reinstall Python from python.org |
| `python: command not found` | Python isn't on your PATH | try `python3`, or reinstall from python.org with "Add to PATH" ticked |

## Requirements

Python 3.8 or newer. Nothing to install — no `pip`, no numpy, no internet.

---

Longer version: `README-core.md` explains what all 20 files do.
Full documentation: `README.md`.
