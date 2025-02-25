---
layout: post
title: "Computing QR Decomposition with Householder Vectors"
date: "2025-02-25 16:06 -0500"
categories: [numerical methods, linear algebra]
tags: python matlab sparse
reading_time: 18
excerpt: How do we compute the QR decomposition of a sparse matrix? I explore
    the differences between two methods of how to compute Householder vectors
    in Python and C++.
---

As I've discussed in a [previous post][sparseqr_blog], 
I have been working through [my own C++/Python implementation][csparse_pp] of
the CSparse routines from Tim Davis's [SuiteSparse][suitesparse] package, as
described in his [book][davis_book]. In Chapter 5 of the book, he describes how
to compute the QR decomposition of a sparse matrix using *Householder
reflections*. As I work through the implementation, I test my code against
output from MATLAB and SciPy. I noticed a difference in output between the code
presented in Davis's book and the code used internally by MATLAB and SciPy. In
this post, I explore the differences between the two methods.


# Householder Reflections
A *Householder reflection* is a linear transformation that reflects a vector
across a plane defined by its normal vector. In this context, we use Householder
reflections to find an orthogonal basis for the column space of a matrix, which
is the basis for the [QR decomposition][qr_wikipedia].

The Householder reflection is defined as

$$
H = I - \beta v v^T
$$

where $I$ is the identity matrix,

$$
\beta = \frac{2}{v^T v}
$$

and $v$ is a vector that defines the plane of reflection. It is convenient to
normalize $v$ such that $v_1 = 1$, or

$$
v = \begin{bmatrix} 1 \\ v_2 \\ \vdots \\ v_N \end{bmatrix},
$$

so that we only need to store the elements $v_2, \ldots, v_N$.

We are interested in finding the Householder vector $v$ that will zero out all
elements below the diagonal of a given column of a matrix. We can then apply
these steps successively to each column to find the QR decomposition of the
matrix, or more directly, the $V$ matrix of Householder vectors, a vector
$\beta$ of scaling factors, and the $R$ matrix of the QR decomposition.

The correct Householder vector $v$ will give a reflection matrix that gives

$$
Hx = \pm \|x\| e_1
$$

when applied to a vector $x$, where $e_1$ is the first column of the identity
matrix. The choice of sign for the reflector is arbitrary, but this step is
where Davis and MATLAB/SciPy differ.


# The Davis Method


# The LAPACK Method


# Numerical Consequences


# Conclusions

For a deeper discussion of numerical analysis, especially as applied to the
Householder version of the QR decomposition, see [Trefethen and
Bau](ref:trefethen_book) and [Golub and Van Loan](ref:GVL). For a more in-depth
discussion of the CSparse routines, see [Davis](ref:davis_book).

### References

<ol>
    <li id="ref:davis_book">Davis, Timothy A. (2006). "Direct Methods for Sparse Linear Systems".</li>
    <li id="ref:trefethen_book">Treftethen, Lloyd and David Bau (1997). "Numerical Linear Algebra".</li>
    <li id="ref:GVL">Golub, Gene H. and Charles F. Van Loan (1996). "Matrix Computations".</li>
</ol>


<!-- ---------------------------------------------------------------------- -->
<!--        Appendix -->
<!-- ---------------------------------------------------------------------- -->
<h2 id="app:A">Appendix: Additional Details on the Derivation</h2>
<div class="app-wrapper">
<label for="app-A" class="app-toggle">Show/Hide Appendix &#8862;</label>
<input type="checkbox" id="app-A" class="app-toggle"/>
<div class="appendix" markdown=1>


</div>  <!-- appendix -->
</div>  <!-- app-wrapper -->
<!-- ---------------------------------------------------------------------- -->
<!--        End Appendix -->
<!-- ---------------------------------------------------------------------- -->


<p class="message" markdown=1>
The entire [source code][python_source] for the figures and algorithms in this
post is available on GitHub.
</p>


<!-- links -->
[suitesparse]: https://github.com/DrTimothyAldenDavis/SuiteSparse
[davis_youtube]: https://www.youtube.com/watch?v=1dGRTOwBkQs
[davis_book]: https://epubs.siam.org/doi/book/10.1137/1.9780898718881
[csparse_pp]: https://github.com/broesler/CppSparse
[qr_wikipedia]: https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections 
[python_source]: https://github.com/broesler/CppSparse/blob/main/python/scripts/householder.py
[sparseqr_blog]: {% post_url 2025-02-25-sparseqr %}
[lapack_dgeqrf]: https://netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html
