# Approximation Techniques Using MATLAB/Octave

## Description
This MATLAB/Octave script implements the **Galerkin method** and other **Weighted Residual Methods** to approximate the solution of differential equations. The method is based on the book:

> *The Infinite Element Method Using MATLAB* by **Young W. Hwon** and **Hyochoong Bang**.

The implemented methods include:
- **Galerkin Method (l = 1)**
- **Least Squares Method (l = 2)**
- **Weighted Residual Method (l = 3)**

The script calculates an approximate function \( u(x) \) as a series expansion and evaluates it at \( x = 0.5 \).

---

## Dependencies
To run this script in **Octave**, install and load the `symbolic` package:

```octave
pkg install -forge symbolic
pkg load symbolic
```

---

## Function Usage
The main function is:

```matlab
u_at_x_05 = aprox_methods(N, l)
```

### Inputs:
- `N`  : Number of terms in the series expansion.
- `l`  : Choice of weighting method.
  - `1` = Galerkin Method
  - `2` = Least Squares Method
  - `3` = Weighted Residual Method

### Outputs:
- `u_at_x_05` : The approximate value of the function at \( x = 0.5 \).

---

## Example Usage
Run the following code to compute an approximation using the Galerkin method with **N = 3** terms:

```matlab
N = 3; % Number of terms in the series
l = 1; % Galerkin method

result = aprox_methods(N, l);
disp(['Result at x = 0.5: ', num2str(result)]);
```

---

## Implementation Details
1. Defines a symbolic function \( u(x) \) as a series expansion:
   \[
   u(x) = \sum_{k=1}^{N} a_k x^k (1 - x)
   \]
2. Computes the residual function:
   \[
   R(x) = \frac{d^2 u}{dx^2} - u + x
   \]
3. Selects a weighting function \( w(x) \) based on the chosen method (`l`).
4. Solves the integral equation:
   \[
   \int_0^1 R(x) w(x) dx = 0
   \]
   to find the unknown coefficients \( a_k \).
5. Substitutes the values of \( a_k \) into \( u(x) \) and evaluates it at \( x = 0.5 \).

---

## References
- Young W. Hwon, Hyochoong Bang. *The Infinite Element Method Using MATLAB*. Wiley, 1999.

