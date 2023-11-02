# Storage for constant values

# Formatting codes
bg_gr = "\x1b[6;30;42m"
bg_bl = "\x1b[0;30;44m"
bg_rd = "\x1b[0;30;41m"
no_bg = "\x1b[0m"
red = "\033[0;31m"
green = "\033[0;32m"
yellow = "\033[0;33m"
blue = "\033[0;34m"
pink = "\033[0;35m"
teal = "\033[0;36m"
gray = "\033[0;37m"
underline = "\033[4m"
bold = "\033[1m"
italic = "\033[3m"
reset = "\033[0m"
line_length = 48

# Formatting values
superscripts = ["\u2070", "\u00B9", "\u00B2", "\u00B3", "\u2074", 
    "\u2075", "\u2076", "\u2077", "\u2078", "\u2079"]   

subscripts = ["\u2080", "\u2081", "\u2082", "\u2083", "\u2084", 
    "\u2085", "\u2086", "\u2087", "\u2088", "\u2089"]  

superscript_dot = "\u02D9"
subscript_dot = "."
superscript_t = "\u1D57"
subscript_t = "\u209C"
superscript_x = "\u02E3"
subscript_x = "\u2093"
superscript_y = "\u02B8"
subscript_y = "\u1D67"
superscript_z = "\u1DBB"
subscript_z = "\u1D69"
superscript_hyp = "\u207B"
subscript_hyp = "\u208B"
half_unicode = "\u00BD"
nabla = "\u2207"

def superscript(term: str):
    result = ""
    for c in term:
        if c >= '0' and c <= '9':
            n = ord(c) - 48
            result += superscripts[n]
        elif c == '.':
            result += superscript_dot
        elif c == "t":
            result += blue + superscript_t + reset
        elif c == "x":
            result += blue + superscript_x + reset
        elif c == "y":
            result += blue + superscript_y + reset
        elif c == "z":
            result += blue + superscript_z + reset
        elif c == "-":
            result += superscript_hyp

    return result

def subscript(term: str):
    result = ""
    for c in term:
        if c >= '0' and c <= '9':
            n = ord(c) - 48
            result += subscripts[n]
        elif c == '.':
            result += subscript_dot
        elif c == "t":
            result += blue + subscript_t + reset
        elif c == "x":
            result += blue + subscript_x + reset
        elif c == "y":
            result += blue + subscript_y + reset
        elif c == "z":
            result += blue + subscript_z + reset
        elif c == "-":
            result += subscript_hyp

    return result

# Linear Algebra symbols
left_br = "\u3008"
right_br = "\u3009"
top_left = "\u2e22"
top_right = "\u2e23"
bottom_left = "\u2e24"
bottom_right = "\u2e25"

# Round floats to how many decimals in printing
decimals = 4
rounding_precision = 0.0000000001

# Trigonmetric functions
TRIG = ["sin", "cos", "tan", "csc", "sec", "cot"]

# Inverse Trigonmetric functions
INV_TRIG = ["arcsin", "arccos", "arctan", "arccsc", "arcsec", "arccot"]

# Numerical constants
EULER = 2.7182818284590452
PI = 3.1415926535897932
pi_str = "\u03c0"
