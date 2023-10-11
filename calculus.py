from constants import *
import copy
import math


"""Abstraction of an algebraic term."""
class Term():
    def __hash__(self):
        return hash(self.__str__())
    
    def __add__(self, other):
        return Sum([self, other])
    
    def __radd__(self, other):
        return Sum([other, self])

    def __mul__(self, other):
        return Product([self, other])

    def derivative(self, show=True):
        c = copy.deepcopy(self)
        d = c.derive(1)
        if show:
            print(f"{italic}f{reset}({blue}x{reset}) = {self}")
            print(f"{italic}f{reset}'({blue}x{reset}) = {d}")
        return d
    
    def integral(self, show=True):
        c = copy.deepcopy(self)
        i = c.integrate(1)
        if show:
            print(f"{italic}f{reset}({blue}x{reset}) = {self}")
            if i is None:
                print(f"{red}(Integral not found){reset}")
            else:
                if str(i)[-12:] != f"{yellow}C{reset}":
                    p1 = f"∫ {italic}f{reset}({blue}x{reset})d"
                    p2 = f"{blue}x{reset} = {i} + {yellow}C{reset}"
                    print(p1 + p2)
                else:
                    print(f"∫ {italic}f{reset}({blue}x{reset})d{blue}x{reset} = {i}")
        return i
    
    def definite_integral(self, a, b, show=True):
        c = copy.deepcopy(self)
        i = c.integral(show=False)
        if isinstance(i, Sum):
            terms = list(filter(lambda t: not isinstance(t, Variable), i.terms))
            i = Sum(terms)

        if show:
            print(f"{italic}f{reset}({blue}x{reset}) = {self}")
            if i is None:
                print(f"{red}(Integral not found){reset}")
                return None
            else:
                sub_a = subscript(str(a))
                sub_b = subscript(str(b))
                sub_term = f"{green}{sub_a}\u2192{sub_b}{reset}"
                p1 = f"∫{sub_term}{italic} f{reset}({blue}x{reset})d"
                p2 = f"{blue}x{reset} = {i}{green}|{reset}{sub_term}"
                print(p1 + p2)

                val_a = i.evaluate(a)
                val_b = i.evaluate(b)
                diff = val_b - val_a
                print(f" = {val_b} - {val_a} = {diff}")
                return diff
        elif i is None:
            return None
        else:
            val_a = i.evaluate(a)
            val_b = i.evaluate(b)
            return val_b - val_a


"""Abstraction of a constant value."""
class Constant(Term):
    def __init__(self, val):
        if (not isinstance(val, int)) and not (isinstance(val, float)):
            raise TypeError("Must be numeric")
        self.val = val        
        self.inner = None

        # Tolerance for grouping to known constants
        self.tolerance = 0.000000001

    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val + other)
        elif isinstance(other, Constant):
            return Constant(self.val + other.val)
        else:
            return Sum([self, other])
    
    def __radd__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val + other)
        elif isinstance(other, Constant):
            return Constant(self.val + other.val)
        else:
            return Sum([other, self])
    
    def __sub__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val - other)
        elif not isinstance(other, Constant):
            print("ERROR: Cannot subtract non-constant")
            return self
        else:
            return Constant(self.val - other.val)
        
    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val * other)
        elif isinstance(other, Constant):
            return Constant(self.val * other.val)
        else:
            return PolyTerm(self.val, 1, inner=other)
    
    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val / other)
        elif isinstance(other, Constant):
            return Constant(self.val / other.val)
        else:
            return PolyTerm(1 / self.val, 1, inner=other)
    
    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val * other)
        elif isinstance(other, Constant):
            return Constant(self.val * other.val)
        else:
            return PolyTerm(self.val, 1, inner=other)

    def __pow__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Constant(self.val ** other)
        elif isinstance(other, Constant):
            return Constant(self.val ** other.val)
        else:
            return ExpTerm(base=self.val, inner=other)

    def __eq__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.val == other
        elif not other.is_constant():
            return False
        else:
            return self.val == other.val

    def __gt__(self, other):
        if not isinstance(other, Constant):
            return False
        else:
            return other.val > self.val
        
    def __lt__(self, other):
        if not isinstance(other, Constant):
            return True
        else:
            return other.val < self.val

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        # Take out zero
        if self.val == 0:
            return "0"

        # Factor out π
        if self.val % PI < self.tolerance:
            coeff = round(self.val / PI)
            return str(coeff) + f"{green}π{reset}" if abs(coeff) > 1 else f"{green}π{reset}"
        
        for denom in range(2, 11):
            if (self.val % (PI / denom)) < self.tolerance:
                coeff = round(self.val / (PI / denom))
                if coeff == 1 and denom == 2:
                    return half_unicode + f" {green}π{reset}"
                return f"{coeff}/{denom} " + f"{green}π{reset}"
        
        # Factor out powers of π        
        for p in range(2, 5):
            if (self.val % (PI ** p)) < self.tolerance:
                pow = round(self.val / (PI ** p))
                exp = superscript(str(pow))
                return f"{green}π{reset}{exp}"

        # Fector out e
        if self.val % EULER < self.tolerance:
            coeff = round(self.val / EULER)
            return str(coeff) + f"{green}e{reset}" if coeff > 1 else f"{green}e{reset}"

        for denom in range(2, 11):
            if (self.val % (EULER / denom)) < self.tolerance:
                coeff = round((self.val / EULER) / denom)
                if coeff == 1 and denom == 2:
                    return half_unicode + f" {green}e{reset}"
                return f"{coeff}/{denom} " + f"{green}e{reset}"

        # Factor out natural logs
        for l in range(2, 11):
            if self.val % math.log(l) < self.tolerance:
                coeff = round(self.val / math.log(l))
                return str(coeff) + f"{green}ln({l}){reset}" if coeff > 1 else f"{green}ln({l}){reset}"

            for denom in range(2, 11):
                if (self.val % (math.log(l) / denom)) < self.tolerance:
                    coeff = round((self.val / math.log(l)) / denom)
                    if coeff == 1 and denom == 2:
                        return half_unicode + f" {green}ln({l}){reset}"
                    return f"{coeff}/{denom} {green}ln({l}){reset}"

        if isinstance(self.val, float):
            return str(round(self.val, decimals))
        return str(self.val)
        
    def flip_sign(self):
        return Constant(-self.val)
    
    def is_constant(self):
        return True
    
    def is_like(self, other):
        return isinstance(other, Constant)
    
    def is_function_of(self, variable):
        return False

    """Takes the nth derivative of term."""
    def derive(self, n):
        return Constant(0)
    
    """Takes the nth partial derivative of term."""
    def partial(self, val, n):
        return Constant(0)

    """Takes the nth integral of term."""
    def integrate(self, n):
        # Base case
        if n == 0:
            return self
        elif self.val == 0:
            return Constant(0)
        elif n == 1:
            return PolyTerm(self.val, 1)
        else:
            return PolyTerm(self.val, 1).integrate(n - 1)
    
    """Takes the nth integral of term with respect to var."""
    def partial_integrate(self, var, n):
        # Base case
        if n == 0:
            return self
        elif self.val == 0:
            return Constant(0)
        elif n == 1:
            return PolyTerm(self.val, 1, variable=var)
        else:
            return PolyTerm(self.val, 1, variable=var).integrate(n - 1)

    def evaluate(self, x):
        return self.val


"""Abstraction of a variable value."""
class Variable(Term):
    def __init__(self, name):
        self.name = name        
        self.inner = None

    def __add__(self, other):
        if not isinstance(other, Variable):
            print("ERROR: Cannot add to non-variable")
            return self
        else:
            return self

    def __eq__(self, other):
        if isinstance(other, Variable):
            return self.name == other.name
        else:
            return False

    def __gt__(self, other):
        if not isinstance(other, Variable) or isinstance(other, Constant):
            return False
        else:
            return True
        
    def __lt__(self, other):
        if not isinstance(other, Variable) or isinstance(other, Constant):
            return True
        else:
            return False

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        return yellow + self.name + reset

    def is_constant(self):
        return False
    
    """Takes the nth derivative of the variable."""
    def derive(self, n):
        return Constant(0)
    
    """Takes the nth integral of the variable."""
    def integrate(self, n):
        # Base case
        if n == 0:
            return self
        else:
            p = PolyTerm(self, 1)
            s = Sum([p, Variable("C")])
            return s.integrate(n - 1)
    
    """Solves algebraically this = opposite."""
    def solve(self, opposite):
        return opposite

    def evaluate(self, x):
        return self.val


"""Abstraction of a product."""
class Product(Term):
    def __init__(self, terms):
        if not isinstance(terms, list):
            raise TypeError("Terms must be list")

        full_terms = []
        for t in terms:
            if isinstance(t, Product):
                full_terms.extend(t.terms())
            else:
                full_terms.append(t)

        self.poly = PolyTerm(1, 0)
        self.other_polys = []
        self.trig = []
        self.logs = []
        self.exps = []
        self.others = []
        for t in full_terms:
            if isinstance(t, Constant):
                self.poly.coeff *= t.val
            
            elif isinstance(t, PolyTerm):
                if t.inner is None:
                    self.poly.coeff *= t.coeff
                    if t.variable == self.poly.variable:
                        self.poly.exp += t.exp
                    else:
                        t.coeff = 1
                        self.other_polys.append(t)

                elif isinstance(t.inner, Product):
                    for t2 in t.inner.terms():
                        if isinstance(t2, PolyTerm):
                            if t2.variable == self.poly.variable:
                                self.poly.exp += (t2.exp * t.exp)
                            else:
                                t2.exp *= t.exp
                                self.other_polys.append(t2)
                        else:
                            new_p = PolyTerm(1, t.exp, inner=t2)
                            self.other.append(new_p)

                elif isinstance(t.inner, TrigTerm):
                    self.poly.coeff *= t.coeff
                    if t.exp == 1:
                        new = t.inner
                    else:
                        new = t
                        new.coeff = 1
                    self.trig.append(new)
                
                elif isinstance(t.inner, ExpTerm):
                    self.poly.coeff *= t.coeff
                    if t.exp == 1:
                        new = t.inner
                    else:
                        new = t
                        new.coeff = 1
                    self.exps.append(new)
            
            elif isinstance(t, TrigTerm):
                found = False
                for i, o in enumerate(self.trig):
                    if isinstance(o, TrigTerm):
                        if (o.fn == t.fn) and (o.inner == t.inner):
                            p = PolyTerm(1, 2, inner=o)
                            self.trig.pop(i)
                            self.trig.append(p)
                            found = True
                    elif isinstance(o, PolyTerm) and isinstance(o.inner, TrigTerm):
                        if (o.inner.fn == t.fn) and (o.inner.inner == t.inner):
                            p = PolyTerm(o.coeff, 1 + o.exp, inner=o.inner)
                            self.trig.pop(i)
                            self.trig.append(p)
                            found = True
                    if t.is_recip(o):
                        self.trig.pop(i)
                
                if not found:
                    self.trig.append(t)
            
            elif isinstance(t, ExpTerm):
                t_copy = copy.deepcopy(t)
                for i, o in enumerate(self.exps):
                    if (o.base == t.base):
                        if isinstance(t_copy.inner, Sum):
                            if o.inner is None:
                                t_copy.inner += PolyTerm(1, 1)
                            else:
                                t_copy.inner += o.inner
                        else:
                            if t_copy.inner is None and o.inner is None:
                                t_copy.inner = PolyTerm(2, 1)
                            elif o.inner is None:
                                t_copy.inner = Sum([t_copy.inner, PolyTerm(1, 1)])
                            else:
                                t_copy.inner = Sum([t_copy.inner, o.inner])
                
                self.exps.append(t_copy)

            elif isinstance(t, LogTerm):
                self.logs.append(t)

        # Clean up other variables
        removes = []
        adds = []
        if len(self.other_polys) > 1:
            for i, o1 in enumerate(self.other_polys[:-1]):
                v1 = o1.variable
                for j, o2 in enumerate(self.other_polys[i + 1:]):
                    if j not in removes:
                        v2 = o2.variable
                        if v1 == v2:
                            new_v = PolyTerm(1, o1.exp + o2.exp, variable=v1)
                            if (o1.exp + o2.exp != 0):
                                adds.append(new_v)
                            removes.append(i)
                            removes.append(i + j + 1)
        
        removes.sort()
        for i in removes[::-1]:
            del self.other_polys[i]
        self.other_polys.extend(adds)

        if self.poly.coeff == 0:
            self.trig = []
            self.logs = []
            self.exps = []
            self.other_polys = []
            self.others = []

    def terms(self):
        if self.poly == PolyTerm(-1, 0) or self.poly == PolyTerm(1, 0):
            t = []
        else:
            t = [self.poly]
        t.extend(self.trig)
        t.extend(self.logs)
        t.extend(self.exps)
        t.extend(self.other_polys)
        t.extend(self.others)
        return t

    def __add__(self, other):
        if other == self:
            return PolyTerm(2, 1, inner=self)
        return Sum([self, other])

    def __sub__(self, other):
        if other == self:
            return Constant(0)
        return Sum([self,PolyTerm(-1, 1, other)])

    def __eq__(self, other):
        if not isinstance(other, Product):
            return False
        if self.poly != other.poly:
            return False
        if set(self.trig) != set(other.trig):
            return False
        if set(self.logs) != set(other.logs):
            return False
        if set(self.exps) != set(other.exps):
            return False
        if set(self.others) != set(other.others):
            return False

        return True

    def __gt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return True
        elif isinstance(other, TrigTerm) or isinstance(other, ExpTerm):
            return True
        elif isinstance(other, LogTerm):
            return True
        elif isinstance(other, PolyTerm):
            return False
        elif isinstance(other, Product):
            if self.poly > other.poly:
                return True
        
        return False

    def __lt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return False
        elif isinstance(other, TrigTerm) or isinstance(other, ExpTerm):
            return False
        elif isinstance(other, LogTerm):
            return False
        elif isinstance(other, PolyTerm):
            return True
        elif isinstance(other, Product):
            if self.poly > other.poly:
                return False
        
        return True

    def __pow__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            c.poly **= other
            for i, term in enumerate(c.trig):
                c.trig[i] = term ** other
            for i, term in enumerate(c.logs):
                c.logs[i] = term ** other
            for i, term in enumerate(c.exps):
                c.exps[i] = term ** other
            for i, term in enumerate(c.others):
                c.others[i] = term ** other
        elif isinstance(other, Constant):
            c.poly **= other.val
            for i, term in enumerate(c.trig):
                c.trig[i] = term ** other.val
            for i, term in enumerate(c.logs):
                c.logs[i] = term ** other.val
            for i, term in enumerate(c.exps):
                c.exps[i] = term ** other.val
            for i, term in enumerate(c.others):
                c.others[i] = term ** other.val
        return c

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        # First, Poly term
        if self.poly.__str__() in ["1", "1.0"]:
            full_str = ""
        elif self.poly.__str__() in ["-1", "-1.0"]:
            full_str = "-"
        else:
            full_str = self.poly.__str__()

        for t in self.other_polys:
            full_str = full_str + t.__str__()

        # Second, Exp terms
        for t in self.exps:
            full_str = full_str + t.__str__()

        # Third, Trig terms
        for t in self.trig:
            full_str = full_str + t.__str__()
        
        # Fourth, Log terms
        for t in self.logs:
            full_str = full_str + t.__str__()
        
        # Then, Other terms
        for t in self.others:
            full_str = full_str + t.__str__()

        return full_str

    def is_constant(self):
        if not self.poly.is_constant():
            return False
        elif self.poly.coeff == 0:
            return True
        for t in self.trig:
            if not t.is_constant():
                return False
        for t in self.logs:
            if not t.is_constant():
                return False
        for t in self.exps:
            if not t.is_constant():
                return False
        for t in self.other_polys:
            if not t.is_constant():
                return False
        for t in self.others:
            if not t.is_constant():
                return False
        
        return True

    """Checks if two terms are reducible."""
    def is_like(self, other):
        if not isinstance(other, Product):
            return False
        if self.poly.exp != other.poly.exp:
            return False
        if set(self.trig) != set(other.trig):
            return False
        if set(self.exps) != set(other.exps):
            return False
        if set(self.others) != set(other.others):
            return False

        return True
    
    def is_function_of(self, variable):
        for term in self.terms():
            if term.is_function_of(variable):
                return True
        return False

    """Adds a term to the product."""
    def multiply(self, other):
        if isinstance(other, Constant):
            self.poly.coeff *= other.val

        elif isinstance(other, PolyTerm):
            if other.inner is None and other.variable == self.poly.variable:
                self.poly.coeff *= other.coeff
                self.poly.exp += other.exp
            
            else:
                for i, t in enumerate(self.trig):
                    found = False
                    if t.inner == other.inner:
                        found = True
                        self.trig[i].exp += other.exp
                        self.trig[i].coeff *= other.coeff

                    if not found:
                        self.trig.append(other)
        
        elif isinstance(other, TrigTerm):
            for i, t in enumerate(self.trig):
                if isinstance(t, TrigTerm) and (t.fn == other.fn):
                    if t.inner == other.inner:
                        p = PolyTerm(1, 2, inner=other)
                        self.trig[i] = p
                        return
                
                elif isinstance(t, PolyTerm) and (t.inner == other):
                    self.trig[i].exp += 1
                    return

            self.trig.append(other)
        
        elif isinstance(other, LogTerm):
            for i, t in enumerate(self.logs):
                if isinstance(t, LogTerm) and (t.base == other.base):
                    if t.inner == other.inner:
                        p = PolyTerm(1, 2, inner=other)
                        self.logs[i] = p
                        return
                
                elif isinstance(t, PolyTerm) and (t.inner == other):
                    self.logs[i].exp += 1
                    return

            self.logs.append(other)

        elif isinstance(other, ExpTerm):
            found = False
            for i, t in enumerate(self.exps):
                if isinstance(t, ExpTerm) and (t.base == other.base):
                    found = True
                    if isinstance(t.inner, Sum):
                        if other.inner is None:
                            t.inner += PolyTerm(1, 1)
                        else:
                            t.inner += other.inner
                    else:
                        if other.inner is None and t.inner is None:
                            t.inner = PolyTerm(2, 1)
                        elif other.inner is None:
                            t.inner = Sum([t.inner, PolyTerm(1, 1)])
                        else:
                            t.inner = Sum([t.inner, other.inner])

            if not found:
                self.exps.append(other)
        
        elif isinstance(other, Product):
            self.multiply(other.poly)
            for t in other.trig:
                self.multiply(t)
            for t in other.logs:
                self.multiply(t)
            for t in other.others:
                self.multiply(t)
    
    """Takes the nth derivative of the product."""
    def derive(self, n):
        # Base case: no derivative
        if n == 0:
            return self
        
        # Product Rule of Differentiation
        terms = []
        full_terms = [self.poly]
        full_terms.extend(self.trig)
        full_terms.extend(self.logs)
        full_terms.extend(self.exps)
        full_terms.extend(self.others)
        for i, term in enumerate(full_terms):
            d = term.derive(1)
            p = d if isinstance(d, Product) else Product([d])
            for j, term_2 in enumerate(full_terms):
                if i != j:
                    p.multiply(term_2)
            terms.append(p)
                
        s = Sum(terms[::-1])
        return s.derive(n - 1)

    """Takes the nth derivative of the product with respect to val."""
    def partial(self, var, n):
        # Base case
        if n == 0:
            return self
        
        relevant_terms = []
        other_terms = []
        for term in self.terms():
            if term.is_function_of(var):
                relevant_terms.append(term)
            else:
                other_terms.append(term)
        
        if len(relevant_terms) == 0:
            return Constant(0)

        if len(relevant_terms) > 1:
            relevant = Product(relevant_terms)
        else:
            relevant = relevant_terms[0]

        relevant_deriv = relevant.partial(var, n)
        other_terms.extend([relevant_deriv])
        if len(other_terms) == 2:
            if other_terms[0] == Constant(1):
                return other_terms[1]
            elif other_terms[1] == Constant(1):
                return other_terms[0]
        return Product(other_terms)

    """Takes the nth integral of the product."""
    def integrate(self, n):
        if n == 0:
            return self
        
        else:
            terms = self.terms()
            if len(terms) == 1:
                return terms[0].integrate(n)

            if len(terms) == 2:
                # Coefficients carry over
                if isinstance(terms[0], PolyTerm) and terms[0].exp == 0:
                    s = terms[1].integrate(n)
                    return s.multiply(terms[0])
                elif isinstance(terms[1], PolyTerm) and terms[1].exp == 0:
                    s = terms[0].integrate(n)
                    return s.multiply(terms[1])
                elif isinstance(terms[0], Constant):
                    s = terms[1].integrate(n)
                    return s.multiply(terms[0])
                elif isinstance(terms[1], Constant):
                    s = terms[0].integrate(n)
                    return s.multiply(terms[1])

                # Integration by Parts
                elif terms[0].inner is None and terms[1].inner is None:
                    # Trig Options
                    if isinstance(terms[0], TrigTerm) and isinstance(terms[1], TrigTerm):
                        if terms[0].fn == "sec" and terms[1].fn == "tan":
                            return Sum([TrigTerm("sec"), Variable("C")])
                        
                        elif terms[1].fn == "sec" and terms[0].fn == "tan":
                            return Sum([TrigTerm("sec"), Variable("C")])
                        
                        elif terms[0].fn == "csc" and terms[1].fn == "cot":
                            return Sum([Product([Constant(-1), TrigTerm("sec")]), Variable("C")])
                        
                        elif terms[1].fn == "csc" and terms[0].fn == "cot":
                            return Sum([Product([Constant(-1), TrigTerm("sec")]), Variable("C")])

                    if isinstance(terms[0], LogTerm):
                        u = terms[0]
                        dv = terms[1]
                    elif isinstance(terms[1], LogTerm):
                        u = terms[1]
                        dv = terms[0]
                    elif isinstance(terms[0], PolyTerm):
                        u = terms[0]
                        dv = terms[1]
                    elif isinstance(terms[1], PolyTerm):
                        u = terms[1]
                        dv = terms[0]
                    elif isinstance(terms[0], ExpTerm):
                        u = terms[0]
                        dv = terms[1]
                    elif isinstance(terms[1], ExpTerm):
                        u = terms[1]
                        dv = terms[0]
                    elif isinstance(terms[0], TrigTerm):
                        u = terms[0]
                        dv = terms[1]
                    elif isinstance(terms[1], TrigTerm):
                        u = terms[1]
                        dv = terms[0]
                    else:
                        # Too complex
                        return None

                    du = copy.deepcopy(u).derive(1)
                    v = copy.deepcopy(dv).integrate(1)
                    v.remove(Variable("C"))
                    if len(v.terms) == 1:
                        v = v.terms[0]

                    uv = Product([u, v])
                    vdu = Product([v, du])
                    vdu.multiply(Constant(-1))
                    return Sum([uv, vdu.integrate(n), Variable("C")])

                elif terms[0].inner is not None and terms[1].inner is None:
                    t0_copy = copy.deepcopy(terms[0])

                    # U Substitution
                    if t0_copy.inner.derive(1).is_like(terms[1]):
                        f = terms[0]
                        g = terms[0].inner
                        g_prime = terms[1]

                        if g_prime == g.derive(1):
                            coeff = None
                        elif isinstance(g_prime, PolyTerm):
                            coeff = g_prime.coeff / g.derive(1).coeff
                        else: 
                            coeff = None
                        
                        f_u = copy.deepcopy(f)
                        f_u.inner = None
                        int_f_u = f_u.integrate(1).multiply(Constant(coeff))
                        i = int_f_u.integrate(n - 1)
                        return Sum([i, Variable("C")])

                    # Integration by parts
                    else:
                        u = terms[1]
                        du = copy.deepcopy(u).derive(1)
                        dv = terms[0]
                        v = copy.deepcopy(dv).integrate(1)
                        uv = Product([u, v])
                        vdu = Product([v, du])
                        vdu.multiply(Constant(-1))
                        return Sum([uv, vdu.integrate(n), Variable("C")])

                elif terms[1].inner is not None and terms[0].inner is None:
                    # U Substitution
                    t1_copy = copy.deepcopy(terms[1])
                    if t1_copy.inner.derive(1).is_like(terms[0]):
                        f = terms[1]
                        g = terms[1].inner
                        g_prime = terms[0]

                        if g_prime == g.derive(1):
                            coeff = None
                        elif isinstance(g_prime, PolyTerm):
                            coeff = g_prime.coeff / g.derive(1).coeff
                        else: 
                            coeff = None
                        
                        f_u = copy.deepcopy(f)
                        f_u.inner = None
                        int_f_u = f_u.integrate(1).multiply(Constant(coeff))
                        i = int_f_u.integrate(n - 1)
                        return Sum([i, Variable("C")])

                    # Integration by parts
                    else:
                        u = terms[0]
                        du = copy.deepcopy(u).derive(1)
                        dv = terms[1]
                        v = copy.deepcopy(dv).integrate(1)
                        uv = u * v
                        vdu = v * du
                        vdu.multiply(Constant(-1))
                        return Sum([uv, vdu.integrate(n), Variable("C")])

                elif isinstance(terms[0], PolyTerm) and isinstance(terms[1], PolyTerm):
                    # Trig Products
                    if isinstance(terms[0].inner, TrigTerm) and isinstance(terms[1].inner, TrigTerm):
                        sincos = False
                        if terms[0].inner.fn == "sin" and terms[1].inner.fn == "cos":
                            sin_term = terms[0]
                            cos_term = terms[1]
                            sincos = True
                        elif terms[0].inner.fn == "cos" and terms[1].inner.fn == "sin":
                            sin_term = terms[1]
                            cos_term = terms[0]
                            sincos = True
                        if sincos:
                            n_ = sin_term.exp
                            m_ = cos_term.exp
                            if n_ == 3:
                                # U Substitution
                                s = Sum([PolyTerm(1, m_), PolyTerm(-1, m_ + 2)])
                                return s.integrate(n)
                            elif m_ == 3:
                                # U Substitution
                                s = Sum([PolyTerm(1, n_), PolyTerm(-1, n_ + 2)])
                                return s.integrate(n)
                            else:
                                # TODO
                                return None
                        
                        tansec = False
                        if terms[0].inner.fn == "tan" and terms[1].inner.fn == "sec":
                            tan_term = terms[0]
                            sec_term = terms[1]
                            tansec = True
                        elif terms[0].inner.fn == "sec" and terms[1].inner.fn == "tan":
                            tan_term = terms[1]
                            sec_term = terms[0]
                            tansec = True
                        if tansec:
                            n_ = tan_term.exp
                            m_ = sec_term.exp
                            if n_ == 3:
                                # U Substitution
                                s = Sum([PolyTerm(1, m_ + 1), PolyTerm(-1, m_ - 1)])
                                return s.integrate(n)
                            elif m_ == 2:
                                # U Substitution
                                s = PolyTerm(1, n_)
                                return s.integrate(n)
                            elif m_ == 4:
                                # U Substitution
                                s = Sum([PolyTerm(1, n_)], PolyTerm(1, 2 + n_))
                                return s.integrate(n)
                            else:
                                # Case by case basis
                                return None

                else:
                    return None

            else:
                return None

    """Takes the nth integral of the product with respect to val."""
    def partial_integrate(self, var, n):
        # Base case
        if n == 0:
            return self
        
        relevant_terms = []
        other_terms = []
        for term in self.terms():
            if term.is_function_of(var):
                relevant_terms.append(term)
            else:
                other_terms.append(term)
        
        if len(relevant_terms) == 0:
            return Product([other_terms, PolyTerm(1, 1, variable=var)])

        if len(relevant_terms) > 1:
            relevant = Product(relevant_terms)
        else:
            relevant = relevant_terms[0]
        relevant_int = relevant.partial_integrate(var, n)

        # Strip constant from sum
        if isinstance(relevant_int, Sum) and (len(relevant_int.terms) == 2):
            relevant_int = relevant_int.terms[0]
        other_terms.extend([relevant_int])
        return Product(other_terms)

    def evaluate(self, x):
        p = 1
        for t in self.terms():
            p *= t.evaluate(x)
        return p


"""Abstraction of a polynomial term in a function."""
class PolyTerm(Term):
    def __init__(self, coeff, exp, inner=None, variable="x"):
        self.coeff = coeff
        self.exp = exp
        self.variable = variable

        # If the inner function is also a polynomial, stack up
        if isinstance(inner, PolyTerm):
            self.coeff *= inner.coeff
            self.exp += inner.exp
            self.inner = inner.inner

        elif isinstance(inner, Product):
            non_constants = []
            for t in inner.terms():
                if isinstance(t, Constant):
                    self.coeff *= inner.coeff ** self.exp
                else:
                    non_constants.append(t)
            self.inner = Product(non_constants)

        else:
            self.inner = inner

        if isinstance(self.coeff, float) and self.coeff.is_integer():
            self.coeff = int(self.coeff)
        if isinstance(self.exp, float) and self.exp.is_integer():
            self.exp = int(self.coeff)

    def __add__(self, other):
        if not isinstance(other, PolyTerm):
            return Sum([self, other])
        
        if self.exp != other.exp:
            # Terms not compatible
            return Sum([self, other])
        
        if self.inner != other.inner:
            # Terms not compatible
            return Sum([self, other])
        
        if self.variable != other.variable:
            # Terms not compatible
            return Sum([self, other])
        
        new_c = self.coeff + other.coeff
        if new_c == 0:
            return Constant(0)

        return PolyTerm(new_c, self.exp, self.inner)
    
    def __sub__(self, other):
        if not isinstance(other, PolyTerm):
            return Sum([self, PolyTerm(-1, 1, inner=other)])
        
        if self.exp != other.exp:
            # Terms not compatible
            return Sum([self, PolyTerm(-1, 1, inner=other)])
        
        if self.inner != other.inner:
            # Terms not compatible
            return Sum([self, PolyTerm(-1, 1, inner=other)])
        
        if self.variable != other.variable:
            # Terms not compatible
            return Sum([self, PolyTerm(-1, 1, inner=other)])
        
        new_c = self.coeff - other.coeff
        if new_c == 0:
            return Constant(0)

        return PolyTerm(new_c, self.exp, self.inner)

    def __pow__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            c.coeff **= other
            if isinstance(c.inner, ExpTerm):
                c.inner **= other
            else:
                c.exp *= other
            return c
        elif isinstance(other, Constant) and other.val.is_integer():
            c.coeff **= other.val
            c.exp *= other.val
            if c.inner is not None:
                c.inner **= other.val
            return c
        else:
            pass

    def __mul__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            c.coeff *= other
            return c
        elif isinstance(other, Constant):
            c.coeff *= other.val
            return c
        elif isinstance(other, PolyTerm):
            if c.inner == other.inner and c.variable == other.variable:
                c.coeff * other.coeff
                c.exp += other.exp
                return c
        return Product([c, other])

    def __truediv__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            c.coeff /= other
        elif isinstance(other, Constant):
            c.coeff /= other.val
        elif isinstance(other, PolyTerm):
            if c.inner == other.inner and c.variable == other.variable:
                c.coeff /= other.coeff
                c.exp -= other.exp
        return c

    def __eq__(self, other):
        if not isinstance(other, PolyTerm):
            return False
        else:
            a = self.coeff == other.coeff
            b = self.exp == other.exp
            c = self.inner == other.inner
            d = self.variable == other.variable
            return all([a, b, c, d])

    def __gt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return True        
        elif isinstance(other, TrigTerm):
            return True
        elif isinstance(other, Product):
            return False
        elif isinstance(other, PolyTerm):
            if self.exp > other.exp:
                return True
            elif self.exp == other.exp:
                return self.coeff > other.coeff
            return False

    def __lt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return False
        elif isinstance(other, TrigTerm):
            return False
        elif isinstance(other, Product):
            return True
        elif isinstance(other, PolyTerm):
            if self.exp > other.exp:
                return False
            elif self.exp == other.exp:
                return self.coeff > other.coeff
            return True
    
    def __str__(self):
        # Base case: 0
        if self.coeff == 0:
            return "0"

        # Base case: constant number
        if self.exp == 0:
            return str(round(self.coeff, decimals))

        ex_str = []
        if isinstance(self.exp, float):
            ex = str(round(self.exp, decimals))
        else:
            ex = str(self.exp)

        for c in ex:
            if c == "-":
                ex_str.append("\u207B")
            elif c == ".":
                ex_str.append(superscript_dot)
            else:
                i = int(c)
                ex_str.append(superscripts[i])

        # Determine coefficient
        if self.coeff == 1:
            c = ""
        elif self.coeff == -1:
            c = "-"
        else:
            if isinstance(self.coeff, float):
                c = str(round(self.coeff, decimals))
            else:
                c = str(self.coeff)

        if self.inner is None:
            if self.exp == 1:
                return f"{c}{blue}{self.variable}{reset}"
            else:
                return f"{c}{blue}{self.variable}{reset}{''.join(ex_str)}"

        elif isinstance(self.inner, TrigTerm):
            if self.exp == 1:
                if self.inner.inner is None:
                    return f"{c}{self.inner.fn}({blue}{self.variable}{reset})"
                else:
                    return f"{c}{self.inner.fn}({self.inner.inner})"
            
            elif self.inner.inner is None:
                return f"{c}{self.inner.fn}{''.join(ex_str)}({blue}{self.variable}{reset})"
            else:
                return f"{c}{self.inner.fn}{''.join(ex_str)}({self.inner.inner})"

        else:
            if self.exp == 1:
                return f"{c}{self.inner.__str__()}"
            else:
                return f"{c}({self.inner.__str__()}){''.join(ex_str)}"
    
    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__str__())

    def is_constant(self):
        if self.inner is not None:
            return self.inner.is_constant()
        return (self.exp == 0) or (self.coeff == 0)
    
    def is_like(self, other):
        if not isinstance(other, PolyTerm):
            return False
        else:
            if self.inner == other.inner:
                return (self.exp == other.exp and self.variable == other.variable)
            else:
                return False

    def is_function_of(self, variable):
        if self.inner is None:
            return self.variable == variable
        else:
            return self.inner.is_function_of(variable)

    def multiply(self, other):
        if isinstance(other, PolyTerm) and self.inner is None and self.variable == other.variable:
            return PolyTerm(self.coeff * other.coeff, self.exp + other.exp)
        
        elif isinstance(other, Constant):
            self.coeff *= other.val
            return self
        
        else:
            return Product([self, other])

    """Takes the nth derivative of the term"""
    def derive(self, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        if self.inner is None:
            if self.is_constant():
                return Constant(0)
            
            elif self.exp == 1:
                if n > 1:
                    return Constant(0)
                else:
                    return Constant(self.coeff)
            else:
                return PolyTerm(self.coeff * self.exp, self.exp - 1, variable=self.variable)

        else:
            # Chain rule
            inner_copy = copy.deepcopy(self.inner)
            g_prime = inner_copy.derive(1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.derive(1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            else:
                f_prime.inner = self.inner

            p = f_prime * g_prime
            return p.derive(n - 1)

    """Takes the nth partial derivative of the term with respect to var"""
    def partial(self, var, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        if self.inner is None:
            if self.is_constant():
                return Constant(0)
            
            if self.variable == var:
                if self.exp == 1:
                    if n > 1:
                        return Constant(0)
                    else:
                        return Constant(self.coeff)
                else:
                    return PolyTerm(self.coeff * self.exp, self.exp - 1, variable = var)
            else:
                return Constant(0)

        else:
            # Chain rule
            inner_copy = copy.deepcopy(self.inner)
            g_prime = inner_copy.partial(var, 1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.partial(var, 1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            else:
                f_prime.inner = self.inner

            p = f_prime * g_prime
            return p.partial(var, n - 1)

    """Takes the nth integral of the term"""
    def integrate(self, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif self.inner is None:
            if self.exp == 0:
                p = PolyTerm(self.coeff, 1)
                s = Sum([p, Variable("C")])
                return s.integrate(n - 1)

            elif self.exp == -1:
                l = LogTerm()
                p = PolyTerm(self.coeff, 1, inner=l)
                s = Sum([p, Variable("C")])
                return s.integrate(n - 1)

            else:
                p = PolyTerm(self.coeff / (self.exp + 1), self.exp + 1)
                s = Sum([p, Variable("C")])
                return s.integrate(n - 1)

        # U-Sub...(?)
        else:
            # Common integrals
            if self.exp == 2 and (self.inner == TrigTerm("sec", inner=None)):
                return TrigTerm("tan").integrate(n - 1)
            elif self.exp == 2 and (self.inner == TrigTerm("csc", inner=None)):
                return Product([Constant(-1), TrigTerm("cot")]).integrate(n - 1)
            elif self.exp == -1 and isinstance(self.inner, Sum):
                if len(self.inner.terms) == 2:
                    t1 = self.inner.terms[0]
                    t2 = self.inner.terms[1]
                    if isinstance(t1, Constant) and (t2 == PolyTerm(1, 1)):
                        return LogTerm(inner=self.inner)
                    elif isinstance(t2, Constant) and (t1 == PolyTerm(1, 1)):
                        return LogTerm(inner=self.inner)

            elif self.exp > 1:
                pass
            else:
                return LogTerm(inner=self.inner).integrate(n - 1)

    """Takes the nth integral of the term"""
    def partial_integrate(self, var, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif (self.inner is None) and (self.variable == var):
            if self.exp == 0:
                p = PolyTerm(self.coeff, 1, variable=var)
                s = Sum([p, Variable("C")])
                return s.partial_integrate(var, n - 1)

            elif self.exp == -1:
                l = LogTerm()
                p = PolyTerm(self.coeff, 1, inner=l, variable=var)
                s = Sum([p, Variable("C")])
                return s.partial_integrate(var, n - 1)

            else:
                p = PolyTerm(self.coeff / (self.exp + 1), self.exp + 1, variable=var)
                s = Sum([p, Variable("C")])
                return s.partial_integrate(var, n - 1)
            
        elif not self.is_function_of(var):
            return Product([self, PolyTerm(1, 1, variable=var)]).partial_integrate(var, n - 1)

        # U-Sub...(?)
        elif (self.variable == var):
            # Common integrals
            if (self.exp == 2) and (self.inner == TrigTerm("sec", inner=None)):
                if (self.inner.variable == var):
                    return TrigTerm("tan").integrate(n - 1)
            elif (self.exp == 2) and (self.inner == TrigTerm("csc", inner=None)):
                return Product([Constant(-1), TrigTerm("cot")]).integrate(n - 1)

            elif self.exp > 1:
                pass
            else:
                return LogTerm(inner=self.inner).integrate(n - 1)

    """Solves algebraically this = opposite."""
    def solve(self, opposite):
        o = opposite / self.coeff
        o = o ** (1 / self.exp)
        if self.inner == None:
            return o
        else:
            return self.inner.solve(o)

    def evaluate(self, x):
        if self.inner is None:
            if isinstance(x, list) or isinstance(x, tuple):
                if self.variable == "x":
                    inner = x[0]
                elif self.variable == "y":
                    inner = x[1]
                elif self.variable == "z":
                    inner = x[2]
                else:
                    inner = x[0]
            else:
                inner = x
        else:
            inner = self.inner.evaluate(x)
        return math.pow(inner, self.exp) * self.coeff


"""Abstraction of an trigonometric term in a function."""
class TrigTerm(Term):
    def __init__(self, fn, inner=None, variable="x"):
        # Confirm correct use
        if fn not in TRIG and fn not in INV_TRIG:
            raise ValueError("Invalid trig function")

        self.fn = fn
        self.inner = inner
        self.variable = variable

    def __add__(self, other):        
        if isinstance(other, TrigTerm):
            if (self.fn == other.fn) and (self.inner == other.inner) and (self.variable == other.variable):
                return PolyTerm(2, 1, inner=self)
            else:
                # Terms not compatible
                return Sum([self, other])
        
        elif isinstance(other, PolyTerm):
            if (other.inner == self) and (other.exp == 1) and (other.inner.variable == self.variable):
                return PolyTerm(other.coeff + 1, 1, self)
        
        else:
            # Terms not compatible
            return Sum([self, other])
    
    def __truediv__(self, other):
        if other == 1:
            return self
        if other == 0:
            return None

        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            return PolyTerm(1 / other, 1, inner=c)
        elif isinstance(other, Constant):
            return PolyTerm(1 / other.val, 1, inner=c)
        elif isinstance(other, TrigTerm) and (other.variable == c.variable):
            if c.fn == other.fn and c.inner == other.inner:
                return Constant(1)
            elif c.fn == "sin" and other.fn == "cos" and c.inner == other.inner:
                return TrigTerm("tan", inner=c.inner)
            elif c.fn == "cos" and other.fn == "sin" and c.inner == other.inner:
                return TrigTerm("cot", inner=c.inner)
        else:
            return Product([self, PolyTerm(1, -1, inner=other)])

    def __pow__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            p = PolyTerm(1, other, inner=c)
            return p
        elif isinstance(other, Constant) and other.val.is_integer():
            p = PolyTerm(1, other.val, inner=c)
            return p
        else:
            pass

    def __eq__(self, other):
        if not isinstance(other, TrigTerm):
            return False
        a = self.fn == other.fn
        b = self.inner == other.inner
        c = self.variable == other.variable
        return all([a, b, c])

    def __gt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return True
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return False
        elif isinstance(other, LogTerm):
            return False
        elif isinstance(other, TrigTerm):
            if self.inner is not None and other.inner is None:
                return True
            return False
    
    def __lt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return False
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return True
        elif isinstance(other, LogTerm):
            return True
        elif isinstance(other, TrigTerm):
            if self.inner is not None and other.inner is None:
                return False
            return True

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        if self.inner is not None:
            return f"{self.fn}({self.inner.__str__()})"
        return f"{self.fn}({blue}{self.variable}{reset})"
    
    def __hash__(self):
        return hash(self.__str__())

    def is_constant(self):
        if self.inner is None:
            return False
        else:
            return self.inner.is_constant()
    
    def is_recip(self, other):
        if self.inner != other.inner or self.variable != other.variable:
            return False
        if (self.fn == "sin") and (other.fn == "csc"):
            return True
        if (self.fn == "csc") and (other.fn == "sin"):
            return True
        if (self.fn == "cos") and (other.fn == "sec"):
            return True
        if (self.fn == "sec") and (other.fn == "cos"):
            return True
        if (self.fn == "tan") and (other.fn == "cot"):
            return True
        if (self.fn == "cot") and (other.fn == "tan"):
            return True
        return False

    def is_like(self, other):
        if not isinstance(other, TrigTerm):
            return False
        else:
            if self.inner == other.inner and self.variable == other.variable:
                return self.fn == other.fn
            else:
                return False

    def is_function_of(self, variable):
        if self.inner is None:
            return self.variable == variable
        else:
            return self.inner.is_function_of(variable)

    """Takes the nth derivative of the term"""
    def derive(self, n):
        # Base case: no derivative
        if n == 0:
            return self

        if self.inner is None:
            if self.fn == "sin":
                t = TrigTerm("cos", inner=None, variable=self.variable)
                return t.derive(n - 1)

            elif self.fn == "cos":
                t = TrigTerm("sin", inner=None, variable=self.variable)
                p = PolyTerm(-1, 1, inner=t)
                return p.derive(n - 1)

            elif self.fn == "tan":
                t = TrigTerm("sec", inner=None, variable=self.variable)
                p = PolyTerm(1, 2, inner=t)
                return p.derive(n - 1)

            elif self.fn == "csc":
                p = PolyTerm(-1, 0)
                t1 = TrigTerm("csc", inner=None, variable=self.variable)
                t2 = TrigTerm("cot", inner=None, variable=self.variable)
                return Product([p, t1, t2])
            
            elif self.fn == "sec":
                t1 = TrigTerm("sec", inner=None, variable=self.variable)
                t2 = TrigTerm("tan", inner=None, variable=self.variable)
                return t1 * t2

            elif self.fn == "cot":
                t = TrigTerm("csc", inner=None, variable=self.variable)
                p = PolyTerm(-1, 2, inner=t)
                return p.derive(n - 1)
        
        # Chain rule
        else:
            g = copy.deepcopy(self.inner)
            g_prime = g.derive(1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.derive(1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            elif isinstance(f_prime, PolyTerm):
                f_prime.inner.inner = self.inner
            else:
                f_prime.inner = self.inner
            p = f_prime * g_prime
            return p.derive(n - 1)

    """Takes the nth partial derivative of the term with respect to var"""
    def partial(self, var, n):
        # Base case: no derivative
        if n == 0:
            return self
        
        if (self.inner is None) and (self.variable == var):
            c = copy.deepcopy(self)
            return c.derive(n)
        
        elif (self.inner is None) and (self.variable != var):
            return Constant(0)
        
        elif (self.inner is not None) and not (self.is_function_of(var)):
            return Constant(0)
        
        else:
            g = copy.deepcopy(self.inner)
            g_prime = g.partial(var, 1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.partial(var, 1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            elif isinstance(f_prime, PolyTerm):
                f_prime.inner.inner = self.inner
            else:
                f_prime.inner = self.inner
            p = f_prime * g_prime
            return p.partial(var, n - 1)

    """Takes the nth integral of the term"""
    def integrate(self, n):
        # Base case: no integral
        if n == 0:
            return self

        if self.inner is None:
            if self.fn == "sin":
                t = TrigTerm("cos", inner=None)
                p = PolyTerm(-1, 1, inner=t)
                s = Sum([p, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "cos":
                t = TrigTerm("sin", inner=None)
                s = Sum([t, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "tan":
                t = TrigTerm("sec", inner=None)
                l = LogTerm(inner=t)
                s = Sum([l, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "csc":
                t1 = TrigTerm("csc", inner=None)
                t2 = TrigTerm("cot", inner=None)
                t3 = PolyTerm(-1, 0, inner=t2)
                s = Sum([t1, t3])
                l = LogTerm(inner=s)
                s2 = Sum([l, Variable("C")])
                return s2.integrate(n - 1)
            
            elif self.fn == "sec":
                t1 = TrigTerm("tan", inner=None)
                t2 = TrigTerm("sec", inner=None)
                s = Sum([t1, t2])
                l = LogTerm(inner=s)
                s2 = Sum([l, Variable("C")])
                return s2.integrate(n - 1)

            elif self.fn == "cot":
                t = TrigTerm("sin", inner=None)
                l = LogTerm(inner=t)
                s = Sum([l, Variable("C")])
                return s.integrate(n - 1)

    """Takes the nth integral of the term with respect to var"""
    def partial_integrate(self, var, n):
        # Base case: no integral
        if n == 0:
            return self

        if (self.inner is None) and (self.variable == var):
            if self.fn == "sin":
                t = TrigTerm("cos", inner=None)
                p = PolyTerm(-1, 1, inner=t)
                s = Sum([p, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "cos":
                t = TrigTerm("sin", inner=None)
                s = Sum([t, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "tan":
                t = TrigTerm("sec", inner=None)
                l = LogTerm(inner=t)
                s = Sum([l, Variable("C")])
                return s.integrate(n - 1)

            elif self.fn == "csc":
                t1 = TrigTerm("csc", inner=None)
                t2 = TrigTerm("cot", inner=None)
                t3 = PolyTerm(-1, 0, inner=t2)
                s = Sum([t1, t3])
                l = LogTerm(inner=s)
                s2 = Sum([l, Variable("C")])
                return s2.integrate(n - 1)
            
            elif self.fn == "sec":
                t1 = TrigTerm("tan", inner=None)
                t2 = TrigTerm("sec", inner=None)
                s = Sum([t1, t2])
                l = LogTerm(inner=s)
                s2 = Sum([l, Variable("C")])
                return s2.integrate(n - 1)

            elif self.fn == "cot":
                t = TrigTerm("sin", inner=None)
                l = LogTerm(inner=t)
                s = Sum([l, Variable("C")])
                return s.integrate(n - 1)
            
        elif (self.inner is None) and (self.variable != var):
            return Product([self, PolyTerm(1, 1, variable=var)])
        
        elif (not self.is_function_of(var)):
            return Product([self, PolyTerm(1, 1, variable=var)])

    """Solves algebraically this = opposite."""
    def solve(self, opposite):
        if self.fn == "sin":
            o = math.asin(opposite)
        elif self.fn == "cos":
            o = math.acos(opposite)
        elif self.fn == "tan":
            o = math.atan(opposite)
        elif self.fn == "csc":
            o = 1 / math.asin(opposite)
        elif self.fn == "sec":
            o = 1 / math.acos(opposite)
        elif self.fn == "cot":
            o = 1 / math.atan(opposite)
    
        elif self.fn == "arcsin":
            o = 1 / math.sin(opposite)
        elif self.fn == "arccos":
            o = math.cos(opposite)
        elif self.fn == "arctan":
            o = math.tan(opposite)
        elif self.fn == "arccsc":
            o = 1 / math.sin(opposite)
        elif self.fn == "arcsec":
            o = 1 / math.cos(opposite)
        elif self.fn == "arccot":
            o = 1 / math.tan(opposite)

        if self.inner == None:
            return o
        else:
            return self.inner.solve(o)

    def evaluate(self, x):
        if self.inner is None:
            if isinstance(x, list) or isinstance(x, tuple):
                if self.variable == "x":
                    inner = x[0]
                elif self.variable == "y":
                    inner = x[1]
                elif self.variable == "z":
                    inner = x[2]
                else:
                    inner = x[0]
            else:
                inner = x
        else:
            inner = self.inner.evaluate(x)

        if self.fn == "sin":
            return math.sin(inner)
        elif self.fn == "cos":
            return math.cos(inner)
        elif self.fn == "tan":
            return math.tan(inner)
        elif self.fn == "csc":
            return 1 / math.sin(inner)
        elif self.fn == "sec":
            return 1 / math.cos(inner)
        elif self.fn == "cot":
            return 1 / math.tan(inner)
    
        elif self.fn == "arcsin":
            return math.asin(inner)
        elif self.fn == "arccos":
            return math.acos(inner)
        elif self.fn == "arctan":
            return math.atan(inner)
        elif self.fn == "arccsc":
            return 1 / math.asin(inner)
        elif self.fn == "arcsec":
            return 1 / math.acos(inner)
        elif self.fn == "arccot":
            return 1 / math.atan(inner)


"""Abstraction of a logarithmic term in a function."""
class LogTerm(Term):
    def __init__(self, base="e", inner=None, variable="x"):
        # Confirm correct use
        if (base != "e") and (base <= 0):
            raise ValueError("Base must be positive number")

        self.base = base
        self.inner = inner
        self.variable = variable
    
    def __add__(self, other):        
        if isinstance(other, LogTerm) and self.variable == other.variable:
            if (self.base == other.base) and (self.inner == other.inner):
                return PolyTerm(2, 1, inner=self)
            else:
                # Terms not compatible
                return Sum([self, other])
        
        elif isinstance(other, PolyTerm):
            if (other.inner == self) and (other.exp == 1):
                return PolyTerm(other.coeff + 1, 1, self)
        
        else:
            # Terms not compatible
            return Sum([self, other])
    
    def __eq__(self, other):
        if not isinstance(other, LogTerm):
            return False
        a = self.base == other.base
        b = self.inner == other.inner
        c = self.variable == other.variable
        return all([a, b, c])
    
    def __gt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return True
        elif isinstance(other, TrigTerm) or isinstance(other, ExpTerm):
            return True
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return False
        elif isinstance(other, LogTerm):
            return self.base < other.base
    
    def __lt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return False
        elif isinstance(other, TrigTerm) or isinstance(other, ExpTerm):
            return False
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return True
        elif isinstance(other, LogTerm):
            return self.base > other.base

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        if self.base == "e":
            if self.inner is not None:
                return f"ln({self.inner.__str__()})"
            return f"ln({blue}{self.variable}{reset})"
        else:
            base_str = []
            b = str(self.base)
            for c in b:
                if c == ".":
                    base_str.append(subscript_dot)
                else:
                    i = int(c)
                    base_str.append(subscripts[i])

            if self.inner is not None:
                return f"log{''.join(base_str)}({self.inner.__str__()})"
            return f"log{''.join(base_str)}({blue}{self.variable}{reset})"

    def __hash__(self):
        return hash(self.__str__())

    def is_constant(self):
        if self.inner is not None:
            return self.inner.is_constant()
        return self.base == 1
    
    def is_like(self, other):
        if not isinstance(other, LogTerm):
            return False
        else:
            if self.inner == other.inner:
                return self.base == other.base
            else:
                return False

    def is_function_of(self, variable):
        if self.inner is None:
            return self.variable == variable
        else:
            return self.inner.is_function_of(variable)

    """Takes the nth derivative of the term"""
    def derive(self, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        if self.inner is None:
            if self.is_constant():
                return Constant(0)
            
            elif self.base == "e":
                p = PolyTerm(1, -1, variable=self.variable)
                return p.derive(n - 1)
            
            else:
                l = LogTerm(inner=Constant(self.base), variable=self.variable)
                p = PolyTerm(1 / l.evaluate(1), -1)
                return p.derive(n - 1)
        
        else:
            # Chain rule
            g_prime = self.inner.derive(1)
            f = self
            f.inner = None
            f_prime = f.derive(1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            else:
                f_prime.inner = self.inner
            
            p = f_prime * g_prime
            return p.derive(n - 1)

    """Takes the nth partial derivative of the term with respect to var"""
    def partial(self, var, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        if (self.inner is None) and (self.variable == var):
            c = copy.deepcopy(self)
            return c.derive(n)
        
        elif (self.inner is None) and (self.variable != var):
            return Constant(0)
        
        elif (self.inner is not None) and not (self.is_function_of(var)):
            return Constant(0)
        
        else:
            g = copy.deepcopy(self.inner)
            g_prime = g.partial(var, 1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.derive(1)
            f_prime.inner = self.inner
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            elif isinstance(f_prime, PolyTerm):
                f_prime.inner.inner = self.inner

            p = f_prime * g_prime
            return p.partial(var, n - 1)

    """Takes the nth integral of the term"""
    def integrate(self, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif self.inner is None:
            p = PolyTerm(1, 1) * LogTerm()
            s = Sum([p, PolyTerm(1, 1), Variable("C")])
            return s

        # U-Sub...(?)
        else:
            pass
    
    """Takes the nth integral of the term with respect to var"""
    def partial_integrate(self, var, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif (self.inner is None) and (self.variable == var):
            p = PolyTerm(1, 1) * LogTerm()
            s = Sum([p, PolyTerm(1, 1), Variable("C")])
            return s
        
        elif (self.inner is None) and (self.variable != var):
            return Product([self, PolyTerm(1, 1, variable=var)])

        # U-Sub...(?)
        else:
            pass
    
    """Solves algebraically this = opposite."""
    def solve(self, opposite):
        if self.base == "e":
            o = EULER ** opposite
        else:
            o = self.base ** opposite

        if self.inner == None:
            return o
        else:
            return self.inner.solve(o)

    def evaluate(self, x):
        if self.inner is None:
            if isinstance(x, list) or isinstance(x, tuple):
                if self.variable == "x":
                    inner = x[0]
                elif self.variable == "y":
                    inner = x[1]
                elif self.variable == "z":
                    inner = x[2]
                else:
                    inner = x[0]
            else:
                inner = x
        else:
            inner = self.inner.evaluate(x)

        if self.base == "e":
            return math.log(inner)
        elif self.base == 1:
            return Constant(0)
        else:
            return math.log(inner, self.base)


"""Abstraction of an exponential term in a function."""
class ExpTerm(Term):
    def __init__(self, base="e", inner=None, variable="x"):
        # Confirm correct use
        if (base != "e") and (base <= 0):
            raise ValueError("Base must be positive number")

        self.base = base
        self.inner = inner
        self.variable = variable

    def __add__(self, other):        
        if isinstance(other, ExpTerm):
            if (self.base == other.base) and (self.inner == other.inner) and (self.variable == other.variable):
                return PolyTerm(2, 1, inner=self)
            else:
                return Sum([self, other])
        
        elif isinstance(other, PolyTerm):
            if (other.inner == self) and (other.exp == 1):
                return PolyTerm(other.coeff + 1, 1, self)
            else:
                return Sum([self, other])

        else:
            return Sum([self, other])
    
    def __radd__(self, other):        
        if isinstance(other, ExpTerm):
            if (self.base == other.base) and (self.inner == other.inner) and (self.variable == other.variable):
                return PolyTerm(2, 1, inner=self)
            else:
                return Sum([other, self])
        
        elif isinstance(other, PolyTerm):
            if (other.inner == self) and (other.exp == 1):
                return PolyTerm(other.coeff + 1, 1, self)
            else:
                return Sum([other, self])

        else:
            return Sum([other, self])

    def __mul__(self, other):
        if other == 0:
            return Constant(0)
        elif other == 1:
            return self

        if isinstance(other, int) or isinstance(other, float):
            return PolyTerm(other, 1, inner=self)
        elif isinstance(other, Constant):
            return PolyTerm(other.val, 1, inner=self)
        
        elif isinstance(other, ExpTerm) and self.variable == other.variable:
            if self.inner is None and other.inner is None:
                return ExpTerm(base=self.base, inner=PolyTerm(2, 1))
            elif self.inner is None:
                i = other.inner + 1
            elif other.inner is None:
                i = self.inner + 1
            else:
                i = self.inner + other.inner
            return ExpTerm(base=self.base, inner=i)
        
        else:
            return Product([self, other])

    def __pow__(self, other):
        c = copy.deepcopy(self)
        if isinstance(other, int) or isinstance(other, float):
            if c.inner is None:
                c.inner = PolyTerm(other, 1)
            else:
                c.inner *= other
        elif isinstance(other, Constant):
            if c.inner is None:
                c.inner = PolyTerm(other.val, 1)
            else:
                c.inner *= other.val
        return c

    def __eq__(self, other):
        if not isinstance(other, ExpTerm):
            return False
        a = self.base == other.base
        b = self.inner == other.inner
        c = self.variable == other.variable
        return all([a, b, c])
    
    def __gt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return True
        elif isinstance(other, TrigTerm):
            return True
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return False
        elif isinstance(other, LogTerm):
            return False
        elif isinstance(other, ExpTerm):
            return self.base < other.base
    
    def __lt__(self, other):
        if isinstance(other, Constant) or isinstance(other, Variable):
            return False
        elif isinstance(other, TrigTerm):
            return False
        elif isinstance(other, PolyTerm) or isinstance(other, Product):
            return True
        elif isinstance(other, LogTerm):
            return True
        elif isinstance(other, ExpTerm):
            return self.base > other.base

    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        b = f"{green}e{reset}" if self.base == "e" else self.base
        if self.inner is not None:
            no_blue = self.inner.__str__().replace(f"{blue}{self.variable}{reset}", self.variable)
            exponent = superscript(no_blue)
            return f"{b}{exponent}"

        exponent = superscript(self.variable)
        return f"{b}{exponent}"

    def __hash__(self):
        return hash(self.__str__())

    def is_constant(self):
        if self.inner is not None:
            return self.inner.is_constant()
        return self.base == 1
    
    def is_like(self, other):
        if not isinstance(other, ExpTerm):
            return False
        else:
            if self.inner == other.inner:
                return self.base == other.base
            else:
                return False

    def is_function_of(self, variable):
        if self.inner is None:
            return self.variable == variable
        else:
            return self.inner.is_function_of(variable)

    """Takes the nth derivative of the term"""
    def derive(self, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        elif self.inner is None:
            if self.is_constant():
                return Constant(0)
            
            elif self.base == "e":
                return self.derive(n - 1)
            
            else:
                # TODO: Make ln(base) more readable
                c = copy.deepcopy(self)
                c *= Constant(math.log(self.base))
                return c.derive(n - 1)
        
        else:
            # Chain rule
            g_prime = self.inner.derive(1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.derive(1)
            if isinstance(f_prime, Sum):
                for i, _ in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            elif isinstance(f_prime, PolyTerm):
                f_prime.inner.inner = self.inner
            else:
                f_prime.inner = self.inner
            p = f_prime * g_prime
            return p.derive(n - 1)

    """Takes the nth partial derivative of the term with respect to var"""
    def partial(self, var, n):
        # Base case: 0th derivative
        if n == 0:
            return self

        if (self.inner is None) and (self.variable == var):
            c = copy.deepcopy(self)
            return c.derive(n)
        
        elif (self.inner is None) and (self.variable != var):
            return Constant(0)
        
        elif (self.inner is not None) and not (self.is_function_of(var)):
            return Constant(0)
        
        else:
            g = copy.deepcopy(self.inner)
            g_prime = g.partial(var, 1)
            f = copy.deepcopy(self)
            f.inner = None
            f_prime = f.partial(var, 1)
            if isinstance(f_prime, Sum):
                for i, t in enumerate(f_prime.terms):
                    f_prime.terms[i].inner = self.inner
            elif isinstance(f_prime, Product):
                f_prime.poly.inner = self.inner
                for i in range(len(f_prime.trig)):
                    f_prime.trig[i].inner = self.inner
                for i in range(len(f_prime.logs)):
                    f_prime.logs[i].inner = self.inner
                for i in range(len(f_prime.others)):
                    f_prime.others[i].inner = self.inner
            elif isinstance(f_prime, PolyTerm):
                f_prime.inner.inner = self.inner
            else:
                f_prime.inner = self.inner
            p = f_prime * g_prime
            return p.partial(var, n - 1)

    """Takes the nth integral of the term"""
    def integrate(self, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif self.inner is None:            
            if self.base == "e":
                return self.integrate(n - 1)
            
            else:
                coeff = 1 / math.log(self.base)
                term = Product([Constant(coeff), self])
                return term.integrate(n - 1)

        # U-Sub...(?)
        else:
            if self.base == "e" and isinstance(self.inner, PolyTerm):
                if self.inner.exp == 1:
                    p = PolyTerm(1 / self.inner.coeff, 1, inner=self)
                    return p.integrate(n - 1)
                else:
                    pass
    
    """Takes the nth partial integral of the term with respect to var"""
    def partial_integrate(self, var, n):
        # Base case: 0th integral
        if n == 0:
            return self

        elif (self.inner is None) and (self.variable == var):         
            if self.base == "e":
                c = copy.deepcopy(self)
                return c.integrate(n - 1)
            
            else:
                coeff = 1 / math.log(self.base)
                term = Product([Constant(coeff), self])
                return term.integrate(n - 1)
            
        elif (self.inner is None) and (self.variable != var):
            return Product([self, PolyTerm(1, 1, variable=var)]).partial_integrate(var, n - 1)

        # U-Sub...(?)
        else:
            if self.base == "e" and isinstance(self.inner, PolyTerm):
                if (self.inner.exp == 1) and (self.inner.is_function_of(var)):
                    p = PolyTerm(1 / self.inner.coeff, 1, inner=self)
                    return p.partial_integrate(var, n - 1)
                elif (self.inner.exp == 1) and (not self.inner.is_function_of(var)):
                    return Product([self, PolyTerm(1, 1, variable=var)]).partial_integrate(var, n - 1)
                else:
                    pass
    
    """Solves algebraically this = opposite."""
    def solve(self, opposite):
        if self.base == "e":
            o = math.log(opposite)
        else:
            o = math.log(opposite, self.base)

        if self.inner == None:
            return o
        else:
            return self.inner.solve(o)

    def evaluate(self, x):
        if self.inner is None:
            if isinstance(x, list) or isinstance(x, tuple):
                if self.variable == "x":
                    inner = x[0]
                elif self.variable == "y":
                    inner = x[1]
                elif self.variable == "z":
                    inner = x[2]
                else:
                    inner = x[0]
            else:
                inner = x
        else:
            inner = self.inner.evaluate(x)

        if self.base == "e":
            return EULER ** inner
        elif self.base == 1:
            return Constant(1)
        else:
            return self.base ** inner


"""Abstraction of a series of linearly connected terms."""
class Sum(Term):
    def __init__(self, terms):
        if not isinstance(terms, list):
            raise TypeError("Terms must be list")
        self.terms = terms
        self.tidy()

    def __contains__(self, other):
        return other in self.terms

    def __eq__(self, other):
        if not isinstance(other, Sum):
            return False
        return set(self.terms) == set(other.terms)

    def __len__(self):
        return len(self.terms)
    
    def __str__(self):
        full_str = ""
        if self.__len__() > 0:
            full_str = self.terms[0].__str__()

            for t in self.terms[1:]:
                full_str += " "
                if isinstance(t, PolyTerm) and t.coeff > 0:
                    full_str += "+ "
                elif isinstance(t, PolyTerm) and isinstance(t.coeff, Variable):
                    full_str += "+ "
                elif isinstance(t, Product) and t.poly.coeff > 0:
                    full_str += "+ "
                elif isinstance(t, TrigTerm) or isinstance(t, LogTerm):
                    full_str += "+ "
                elif isinstance(t, ExpTerm) or isinstance(t, Variable):
                    full_str += "+ "
                elif isinstance(t, Constant) and t.val > 0:
                    full_str += "+ "
                full_str += t.__str__()
        else:
            full_str = "0"

        return full_str
    
    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if isinstance(other, Sum):
            self.terms.extend(other.terms)
        
        else:
            self.terms.append(other)
        self.tidy()
        return self

    def __radd__(self, other):
        if isinstance(other, Sum):
            self.terms.extend(other.terms)
        
        else:
            self.terms.append(other)
        self.tidy()
        return self
    
    def __sub__(self, other):
        if other in self.terms:
            self.terms.remove(other)
            return self

        if isinstance(other, Sum):
            for term in other.terms:
                self.terms.append(PolyTerm(-1, 1, inner=term))
        
        else:
            self.terms.append(PolyTerm(-1, 1, inner=other))
        self.tidy()
        return self

    def __pow__(self, other):
        if self.is_constant():
            v = self.evaluate(1)
            if isinstance(other, int) or isinstance(other, float):
                return Constant(v ** other)
            elif isinstance(other, Constant):
                return Constant(v ** other.val)
        else:
            if isinstance(other, int) or isinstance(other, float):
                return PolyTerm(1, other, inner=self)
            elif isinstance(other, Constant):
                return PolyTerm(1, other.val, inner=self)

    def is_constant(self):
        return all([t.is_constant() for t in self.terms])

    def is_function_of(self, variable):
        for term in self.terms:
            if term.is_function_of(variable):
                return True
        return False

    def multiply(self, other):
        if isinstance(other, Sum):
            # Too complicated for now...
            return None
        
        else:
            new_terms = []
            for t in self.terms:
                if isinstance(t, Variable):
                    new_terms.append(t)
                else:
                    p = other * t
                    new_terms.append(p)
            self.terms = new_terms
            self.tidy()
        
        return self

    def remove(self, other):
        if other in self.terms:
            self.terms.remove(other)
        self.tidy()

    """Takes the nth derivative of the sum"""
    def derive(self, n):
        # Base case: no derivative
        if n == 0:
            return self

        terms = self.terms
        for i in range(n):
            terms = [t.derive(1) for t in terms]

        s = Sum(terms)
        s.tidy()
        if len(s.terms) == 1:
            return s.terms[0]
        return s
    
    """Takes the nth partial derivative of the sum wtih respect to var"""
    def partial(self, var, n):
        # Base case: no derivative
        if n == 0:
            return self

        terms = self.terms
        for i in range(n):
            terms = [t.partial(var, 1) for t in terms]
        s = Sum(terms)
        s.tidy()
        return s
    
    """Takes the nth integral of the sum"""
    def integrate(self, n):
        # Base case: no integral
        if n == 0:
            return self

        terms = self.terms
        for i in range(n):
            terms = [t.integrate(1) for t in terms]
        terms.append(Variable("C"))
        s = Sum(terms)
        s.tidy()
        return s
    
    """Takes the integral of the sum with respect to var"""
    def partial_integrate(self, var, n):
        # Base case: no integral
        if n == 0:
            return self

        terms = self.terms
        terms = [t.partial_integrate(var, n) for t in terms]
        terms.append(Variable("C"))
        s = Sum(terms)
        s.tidy()
        return s

    def evaluate(self, x):
        total = 0
        for t in self.terms:
            total += t.evaluate(x)
        return total

    def tidy(self):
        # Clean up sums
        clean_ups = []
        for i, t in enumerate(self.terms):
            if isinstance(t, Sum):
                clean_ups.append(i)
                self.terms.extend(t.terms)
        
        for c in clean_ups[::-1]:
            self.terms.pop(c)

        # Look for Trig identities
        for i, t1 in enumerate(self.terms):
            for j, t2 in enumerate(self.terms):
                if isinstance(t1, PolyTerm) and isinstance(t2, PolyTerm):
                    if isinstance(t1.inner, TrigTerm) and isinstance(t2.inner, TrigTerm):
                        if (t1.exp == 2) and (t2.exp == 2):
                            if set([t1.inner.fn, t2.inner.fn]) == set(["sin", "cos"]):
                                pairs = min([t1.coeff, t2.coeff])
                                remaining = max([t1.coeff, t2.coeff]) - pairs
                                if t1.coeff == t2.coeff:
                                    self.terms[i] = Constant(pairs)
                                    self.terms.pop(j)
                                elif t1.coeff > t2.coeff:
                                    self.terms[j] = Constant(pairs)
                                    self.terms[i].coeff = remaining
                                elif t2.coeff > t1.coeff:
                                    self.terms[i] = Constant(pairs)
                                    self.terms[j].coeff = remaining

        # Clean up constants
        const_ind = -1
        clean_ups = []
        for i, t in enumerate(self.terms):
            if t is None:
                clean_ups.append(i)

            elif t.is_constant():
                if t.evaluate(1) == 0:
                    clean_ups.append(i)
                elif const_ind == -1:
                    const_ind = i
                else:
                    self.terms[const_ind] += t.val
                    clean_ups.append(i)

        for c in clean_ups[::-1]:
            self.terms.pop(c)

        # Clean up variables
        seen = False
        clean_ups = []
        for i, t in enumerate(self.terms):
            if isinstance(t, Variable) and not seen:
                seen = True
            elif isinstance(t, Variable):
                clean_ups.append(i)

        for c in clean_ups[::-1]:
            self.terms.pop(c)

        # Combine like terms
        for i, t1 in enumerate(self.terms):
            for j, t2 in enumerate(self.terms):
                if isinstance(t1, Product) and isinstance(t2, Product):
                    if t1.is_like(t2) and (j > i):
                        self.terms[i].poly.coeff += t2.poly.coeff
                        self.terms.pop(j)

        # Sort by degree
        self.terms.sort(reverse=True)
