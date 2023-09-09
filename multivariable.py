from __future__ import annotations
from calculus import *
from constants import *
import copy
import math


class Vector():
    def __init__(self, vals: list):
        self.vals = vals
        for i, v in enumerate(vals):
            if isinstance(v, int) or isinstance(v, float):
                self.vals[i] = Constant(v)

        self.dimensions = len(vals)
        if self.dimensions == 3:
            self.i = self.vals[0]
            self.j = self.vals[0]
            self.k = self.vals[0]
 
    def __add__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] += other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] += other.val
            return results
        elif isinstance(other, Vector):
            if other.dimensions != self.dimensions:
                raise ValueError("Dimension mismatch in vector addition.")
            else:
                for i in range(self.dimensions):
                    results.vals[i] += other.vals[i]
                return results

    def __sub__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] -= other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] -= other.val
            return results
        elif isinstance(other, Vector):
            if other.dimensions != self.dimensions:
                raise ValueError("Dimension mismatch in vector subtraction.")
            else:
                for i in range(self.dimensions):
                    results.vals[i] -= other.vals[i]
                return results

    def __mul__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] *= other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] *= other.val
            return results
        else:
            pass
    
    def __truediv__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] /= other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] /= other.val
            return results
        else:
            pass

    def __floordiv__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] //= other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] //= other.val
            return results
        else:
            pass

    def __mod__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] %= other
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] %= other.val
        return results
    
    def __pow__(self, other):
        results = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            for i in range(self.dimensions):
                results.vals[i] **= other
            return results
        elif isinstance(other, Constant):
            for i in range(self.dimensions):
                results.vals[i] **= other.val
            return results
        else:
            pass

    def __repr__(self):
        return self.__str__()
    
    def __str__(self, type="standard"):
        if type == "standard":
            return self.standard_row()
        elif type == "column":
            return self.column_vector()
        elif type == "row":
            return self.row_vector()
        elif type == "basis":
            return self.basis_vectors()

    def contents_text(self):
        vals_text = []
        for val in self.vals:
            if isinstance(val, float):
                vals_text.append(str(round(val, decimals)))
            elif isinstance(val, int) or isinstance(val, Constant):
                vals_text.append(str(val))
            else:
                vals_text.append(val.__str__())
        return vals_text

    def standard_row(self):
        result = left_br
        vals_text = self.contents_text()
        result += ", ".join(vals_text)
        result += right_br
        return result
    
    def column_vector(self):
        result = "["
        vals_text = self.contents_text()
        result += "\n ".join(vals_text)
        result += "]"
        return result
    
    def row_vector(self):
        result = "["
        vals_text = self.contents_text()
        result += "  ".join(vals_text)
        result += "]"
        return result
    
    def basis_vectors(self):
        vals_text = self.contents_text()
        if self.dimensions == 2:
            t1 = vals_text[0] + bold + teal + "i " + reset
            t2 = vals_text[1] + bold + teal + "j" + reset
            if self.vals[1].is_constant() and self.vals[1] < 0:
                return f"{t1} {t2}"
            return f"{t1} + {t2}"

        elif self.dimensions == 3:
            t1 = vals_text[0] + bold + teal + "i " + reset
            t2 = vals_text[1] + bold + teal + "j " + reset
            t3 = vals_text[2] + bold + teal + "k" + reset
            result = t1
            if self.vals[1].is_constant() and self.vals[1] < 0:
                result += f" {t2}"
            else:
                result += f"+ {t2}"
            
            if self.vals[2].is_constant() and self.vals[2] < 0:
                result += f" {t3}"
            else:
                result += f"+ {t3}"
            return result
        else:
            pass

    """Returns the magnitude of the vector."""
    def magnitude(self) -> float:
        terms = [v ** 2 for v in self.vals]
        s = Sum(terms)
        return s ** 0.5

    """Returns the unit vector (scaled to size 1) of the vector direciton."""
    def unit_vector(self) -> float:
        mag = self.magnitude()
        if mag is None:
            return None
        return self / mag

    """Returns the dot product of this vector and other."""
    def dot_product(self, other: Vector) -> float:
        if not isinstance(other, Vector):
            return TypeError("Cannot take dot product without two vectors.")
        result = copy.deepcopy(self)
        for i in range(self.dimensions):
            result.vals[i] *= other.vals[i]

        vals = []
        for val in result.vals:
            if isinstance(val, int) or isinstance(val, float):
                vals.append(Constant(val))
            else:
                vals.append(val)

        return sum(vals)

    """Returns the angle theta between this vector and other"""
    def theta(self, other: Vector) -> float:
        numerator = self.dot_product(other)
        denominator = self.magnitude() * other.magnitude()
        if denominator is None or denominator == 0:
            return None
        return math.acos(numerator / denominator)

    """Returns the scalar projection of this vector onto other."""
    def scalar_proj(self, other: Vector) -> float:
        if not isinstance(other, Vector):
            return TypeError("Cannot project onto scalar.")
        if other.dimensions != self.dimensions:
            raise ValueError("Dimension mismatch in projection.")
        num = self.dot_product(other)
        return num / other.magnitude()

    """Returns the vector projection of this vector onto other."""
    def vector_proj(self, other: Vector) -> Vector:
        if not isinstance(other, Vector):
            return TypeError("Cannot project onto scalar.")
        if other.dimensions != self.dimensions:
            raise ValueError("Dimension mismatch in projection.")
        left = self.dot_product(other) / other.magnitude()
        right = other.unit_vector()
        return right * left
    
    """Returns the cross product of this vector and other."""
    def cross_product(self, other: Vector) -> Vector:
        if other.dimensions != self.dimensions:
            raise ValueError("Dimension mismatch in cross product.")
        if self.dimensions == 3:
            t1 = self.vals[1] * other.vals[2]
            t1 -= self.vals[2] * other.vals[1]
            t2 = self.vals[2] * other.vals[0]
            t2 -= self.vals[0] * other.vals[2]
            t3 = self.vals[0] * other.vals[1]
            t3 -= self.vals[1] * other.vals[0]
            return Vector([t1, t2, t3])
        elif self.dimensions == 2:
            t1 = self.vals[0] * other.vals[1]
            t1 -= self.vals[1] * other.vals[0]
            return Vector([0, 0, t1])
        else:
            raise ValueError("Cross Product only 3-Dimensional.")

    """Takes the derivative of a one-variable vector."""
    def derive(self, n=1) -> Vector:
        for i, v in enumerate(self.vals):
            if isinstance(v, int) or isinstance(v, float):
                d = 0
            else:
                d = v.derive(n)
            self.vals[i] = d
        return self
    
    def derivative(self, show=True) -> Vector:
        c = copy.deepcopy(self)
        for i, v in enumerate(c.vals):
            if isinstance(v, int) or isinstance(v, float):
                d = 0
            else:
                d = v.derive(1)
            c.vals[i] = d

        if show:
            print(f"{italic}f{reset}({blue}x{reset}) = {self}")
            print(f"{italic}f{reset}'({blue}x{reset}) = {c}")
        return c
    
    """Takes the integral of a one-variable vector."""
    def integrate(self, n=1) -> Vector:
        for i, v in enumerate(self.vals):
            if isinstance(v, int) or isinstance(v, float):
                t = Constant(v)
                d = t.integrate(n)
            else:
                d = v.integrate(n)
            self.vals[i] = d
        return self

    def integral(self, show=True) -> Vector:
        c = copy.deepcopy(self)
        for i, v in enumerate(c.vals):
            if isinstance(v, int) or isinstance(v, float):
                const = Constant(v)
                d = const.integrate(1)
            else:
                d = v.integrate(1)

            if str(d)[-12:] != f"{yellow}C{reset}":
                c.vals[i] = Sum([d, Variable("C")])
            else:
                c.vals[i] = d

        if show:
            print(f"{italic}f{reset}({blue}x{reset}) = {self}")
            print(f"∫ {italic}f{reset}({blue}x{reset})d{blue}x{reset} = {c}")
        return c

    def unit_tangent_vector(self):
        d = self.derivative(show=False)
        m = d.magnitude()
        T = d / m
        return T

    def unit_normal_vector(self):
        T_prime = self.unit_tangent_vector().derivative(show=False)
        m = T_prime.magnitude()
        return T_prime / m

    def unit_binormal_vector(self):
        T = self.unit_tangent_vector()
        N = self.unit_normal_vector()
        return T.cross_product(N)

    """Calculates the length of the plane or space curve described by r(t)."""
    def curve_length(self, a, b) -> float:
        terms = []
        for v in self.vals:
            c = copy.deepcopy(v)
            t = c.derive(1)
            terms.append(t ** 2)
        s = Sum(terms)
        p = s ** 0.5

        return p.definite_integral(a, b, show=False)

    def curvature(self):
        d = self.derivative(show=False)
        m = d.magnitude()
        T = d / m
        T_prime = T.derivative(show=False)
        return T_prime.magnitude() / m

    def evaluate(self, x):
        terms = []
        for val in self.vals:
            if isinstance(val, int) or isinstance(val, float):
                terms.append(val)
            else:
                terms.append(val.evaluate(x))
        return Vector(terms)


def Gradient(function, dimensions=3) -> Vector:
    if isinstance(function, Vector):
        raise TypeError("Can only take gradient of scalar function.")
    
    if dimensions == 2:
        c1 = copy.deepcopy(function)
        dx = c1.partial("x", 1)
        c2 = copy.deepcopy(function)
        dy = c2.partial("y", 1)
        v = Vector([dx, dy])
        print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {function}")
        print(f"{nabla} {italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {v}")
        return v
    
    elif dimensions == 3:
        c1 = copy.deepcopy(function)
        dx = c1.partial("x", 1)
        c2 = copy.deepcopy(function)
        dy = c2.partial("y", 1)
        c3 = copy.deepcopy(function)
        dz = c3.partial("z", 1)
        v = Vector([dx, dy, dz])
        if function.is_function_of("z"):
            print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {function}")
            print(f"{nabla} {italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {v}")
        else:
            print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {function}")
            print(f"{nabla} {italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {v}")
        return v

def DirectionalDerivative(function, vector: Vector, point, dimensions=3):
    grad = Gradient(function, dimensions=dimensions)
    grad_point = grad.evaluate(point)
    print(f" @ {point}: {grad_point}")
    unit = vector.unit_vector()
    print(f"\n{bold}u{reset} = {unit}")

    dd = grad_point.dot_product(unit)
    print(f"{nabla} {italic}f{reset}{point} • {bold}u{reset} = {dd}")
    return dd

def DoubleInterval(function, limits=None):
    i_x = function.partial_integrate("x", 1)
    if limits is not None:
        x_min = limits[0][0]
        x_max = limits[0][1]
        i_x = i_x.evaluate(x_max, "x") - i_x.evaluate(x_min, "x")
    i_y = i_x.partial_integrate("y", 1)
    if limits is not None:
        y_min = limits[1][0]
        y_max = limits[1][1]
        i_y = i_y.evaluate(y_max, "y") - i_y.evaluate(y_min, "y")
    return i_y

def Divergence(function: Vector, dimensions=3):
    if dimensions == 2:
        f1 = copy.deepcopy(function)
        dx = f1.vals[0].partial("x", 1)
        f2 = copy.deepcopy(function)
        dy = f2.vals[1].partial("y", 1)
        s = dx + dy
        print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {function}")
        print(f"{nabla}⋅{italic}f{reset}({blue}x{reset}, {blue}y{reset}) = {s}")
        return s

    if dimensions == 3:
        f = copy.deepcopy(function)
        dx = f.vals[0].partial("x", 1)
        dy = f.vals[1].partial("y", 1)
        dz = f.vals[2].partial("z", 1)
        s = dx + dy + dz
        print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {function}")
        print(f"{nabla}⋅{italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {s}")
        return s

def Curl(function: Vector) -> Vector:
    f1 = copy.deepcopy(function)
    Ry = f1.vals[2].partial("y", 1)
    Qz = f1.vals[1].partial("z", 1)
    print(f"Ry: {Ry}, Qz: {Qz}")

    f2 = copy.deepcopy(function)
    Pz = f2.vals[0].partial("z", 1)
    Rx = f2.vals[2].partial("x", 1)

    f3 = copy.deepcopy(function)
    Qx = f3.vals[1].partial("x", 1)
    Py = f3.vals[0].partial("y", 1)
    v = Vector([Ry - Qz, Pz - Rx, Qx - Py])
    print(f"{italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {function}")
    print(f"{nabla}×{italic}f{reset}({blue}x{reset}, {blue}y{reset}, {blue}z{reset}) = {v}")
    return v


if __name__ == "__main__":
    e1 = ExpTerm()
    trig1 = TrigTerm("sin", variable="y")
    log1 = LogTerm(variable="z")
    s2 = Sum([e1, trig1, log1])
    #Gradient(s2)

    t1 = PolyTerm(1, 1, variable="x")
    t2 = Product([PolyTerm(1, 1, variable="y"), PolyTerm(1, 1, variable="z")])
    t3 = LogTerm(inner=t2)
    p = Product([t1, t3])
    #Gradient(p)
    
    p1 = Product([PolyTerm(1, 1), PolyTerm(1, 1, variable="z")])
    p2 = Product([PolyTerm(1, 1), PolyTerm(1, 1, variable="y"), PolyTerm(1, 1, variable="z")])
    t3 = PolyTerm(-1, 2, variable="y")
    v = Vector([p1, p2, t3])
    Divergence(v)
    v = Vector([p1, p2, t3])
    Curl(v)
