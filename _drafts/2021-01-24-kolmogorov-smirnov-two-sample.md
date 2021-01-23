---
layout: post
title:  "Kolmogorov-Smirnov Test for Two Samples"
date:   2021-01-23 00:26 -0500
categories: statistics
tags: hypothesis-testing, statistics, python
---
Consider two independent samples $X_1, \dots, X_n$, and
$Y_1, \dots, Y_m$ of independent, real-valued, continuous random
variables, and assume that the $X_i$'s are i.i.d. with some cdf $F$ and
that the $Y_i$'s are i.i.d. with some cdf $G$. Note that the two samples
may have different sizes (if $n \ne m$). We want to test whether
$F = G$. Consider the following hypotheses:

$$
\begin{aligned}
  H_0 \colon ``F = G" \\
  H_1 \colon ``F \ne G"\end{aligned}
$$

For simplicity, we will assume
that $F$ and $G$ are continuous and increasing.

## Example Experiment

An example experiment in which testing if two samples are from the same
distribution is of interest may be encountered in a lab setting where we
have two devices for measurement, and wish to determine if the errors
have the same distribution for our analysis.

## CDF Distributions

Let $$\begin{aligned}
  U_i &= F(X_i), \quad \forall i = 1, \dots, n, \\
  V_j &= G(Y_j), \quad \forall j = 1, \dots, n.\end{aligned}$$

::: {.prop}
**Proposition 1**. *The distribution of the cdf of a continuous random
variable is uniform on $[0,
  1]$.*
:::

::: {.proof}
*Proof.* The distributions of $U_i$ and $V_j$ can be determined by
finding their cdfs. The cdf of $U_i$ is defined by
$F_U(t) \coloneqq \mathbb{P}\left[U_i \le t\right]$. Assuming that
$F(X)$ and $G(Y)$ are invertible, it follows that $$\begin{aligned}
{3}
  \mathbb{P}\left[U_i \le t\right] &= \mathbb{P}\left[F(X_i) \le t\right] &\quad&\text{(definition of $U_i$)} \\
                   &= \mathbb{P}\left[X_i \le F^{-1}(t)\right] \\
                   &= F(F^{-1}(t)) &\quad&\text{(definition of cdf)} \\
                   &= t \\
  \therefore F_U(t) &= t \\
  \implies f_U(t) &= \mathcal{U}\left(\left[ 0, 1 \right]\right) \tag*{\qedhere}\end{aligned}$$ ◻
:::

Likewise, $f_V(t) = \mathcal{U}\left(\left[ 0, 1 \right]\right)$.

## Empirical CDFs

Let $F_n$ be the empirical cdf of $\{X_1, \dots, X_n\}$ and $G_m$ be the
empirical cdf of $\{Y_1, \dots, Y_m\}$.

### The Test Statistic

Let $$T_{n,m} = \sup_{t \in \mathbb{R}} \left| F_n(t) - G_m(t) \right|$$

::: {.prop}
**Proposition 2**. *The test statistic $T_{n,m}$ can be written as the
maximum value of a finite set of numbers.*
:::

::: {.proof}
*Proof.* By definition, the cdf $$\begin{aligned}
{3}
    F(t) &= \mathbb{P}\left[X \le t\right] \quad \forall t \in \mathbb{R}\\
         &= \mathbb{E}\left[\mathbbm{1}\!\left\{X \le t\right\}\right]. \\
    \intertext{By the Law of Large Numbers, the expectation can be approximated
        by the sample average, so we can define the \emph{empirical cdf} as}
    F_n(t) &= \frac{1}{n}\sum_{i=1}^{n} \mathbbm{1}\!\left\{X_i \le t\right\} \addtocounter{equation}{1}\tag{\theequation}\label{eq:F_n}
    \intertext{Likewise,}
    G_m(t) &= \frac{1}{m}\sum_{j=1}^{m} \mathbbm{1}\!\left\{Y_j \le t\right\}. \addtocounter{equation}{1}\tag{\theequation}\label{eq:G_m}
  \end{aligned}$$
$$\therefore T_{n,m} = \sup_{t \in \mathbb{R}} \left| \frac{1}{n}\sum_{i=1}^{n} \mathbbm{1}\!\left\{X_i \le t\right\} - \frac{1}{m}\sum_{j=1}^{m} \mathbbm{1}\!\left\{Y_j \le t\right\} \right|.$$
The empirical cdfs [\[eq:F_n\]](#eq:F_n){reference-type="eqref"
reference="eq:F_n"} and [\[eq:G_m\]](#eq:G_m){reference-type="eqref"
reference="eq:G_m"} can also be written $$\begin{aligned}
{3}
    F_n(t) &= \#\{i=1, \dots, n \colon X_i \le t\} \cdot \frac{1}{n} \\
    G_m(t) &= \#\{i=1, \dots, m \colon Y_j \le t\} \cdot \frac{1}{m},
  \end{aligned}$$ so the only values that the empirical cdfs can take
are the discrete sets $$\begin{aligned}
    F_n(i) &= \frac{i}{n} \quad \forall i = 1, \dots, n \\
    G_m(j) &= \frac{j}{m} \quad \forall j = 1, \dots, m.
  \end{aligned}$$ Therefore, the test statistic can be rewritten as the
maximum value of a finite set of numbers: $$\begin{split}
      T_{n,m} = \max_{i=0,\dots,n} \Bigg[
      &\max_{j=0,\dots,m} \left| \frac{i}{n} - \frac{j}{m} \right| 
        \mathbbm{1}\!\left\{Y^{(j)} \le X^{(i)} < Y^{(j+1)}\right\}, \\ 
      &\max_{k=j+1, \dots, m} \left| \frac{i}{n} - \frac{k}{m} \right| 
        \mathbbm{1}\!\left\{Y^{(k)} \le X^{(i+1)}\right\} \Bigg]
    \end{split}$$ where $X^{(i)}$ is the $i^\text{th}$ value in the
ordered set of data $X^{(1)} \le \cdots \le X^{(n)}$. The values
$X^{(0)}, Y^{(0)} \coloneqq -\infty$ are prepended to the otherwise
finite realizations to simplify the computation. ◻
:::

The following algorithm calculates the KS test statistic for two given
samples.

$X, Y$ are vectors of real numbers. $0 \le T_{n,m} \le 1$.
$X_s \gets \{-\infty,$ $\}$ $Y_s \gets$ $n \gets \dim X_s$
$m \gets \dim Y_s$ $T_v \gets$ empty array of size $n$ $j \gets j$ +
$k \gets j$ + $\displaystyle{T_v^{(i)} \gets 
        \max\left(\left|\frac{i}{n} - \frac{j}{m}\right|,
              \left|\frac{i}{n} - \frac{k}{m}\right|\right)}$
$\max_i T_v$ $A$ is sorted in ascending order.
$\#\{i=1,\dots,\dim A \colon k < A_i\}$

The following subroutine is an implementation of
Algorithm [\[alg:ks_stat\]](#alg:ks_stat){reference-type="ref"
reference="alg:ks_stat"}. It computes an array of values $T_v(i)$ for
each value of $X_i$. The test statistic $T_{n,m}$ is the maximum of
these values.

``` {.python language="python" rangeprefix="\\#\\ <<" rangesuffix=">>" includerangemarker="false" linerange="begin__ks_2samp-end__ks_2samp"}
def _ks_2samp(X, Y):
    """Compute the Kolmogorov-Smirnov statistic on 2 samples.

    Parameters
    ----------
    X, Y : (N,), (M,) array_like
        Two arrays of sample observations assumed to be drawn from a continuous
        distribution, sample sizes can differ.

    Returns
    -------
    Tv : (N+1,) ndarray
        Maximum difference in CDFs for each value of X.

    .. note:: The KS statistic itself is the maximum of these `Tv` values, but
        use this helper function for debugging.
    """
    n = len(X)
    m = len(Y)
    # Sort copies of the data
    Xs = np.hstack([-np.inf, np.sort(X)])  # pad extra point
    Ys = np.sort(Y)
    # Calculate the maximum difference in the empirical CDFs
    Tv = np.zeros(n+1)  # extra value at Fn = 0 (Xs -> -infty)
    js = np.zeros(n+1, dtype=int)
    j = 0
    for i in range(n+1):
        # Find greatest Ys point j s.t. Ys[j] <= Xs[i] and Xs[i] < Ys[j+1]
        j += _rank(Ys[j:], Xs[i])  # only search remaining values

        test_lo = np.abs(i/n - j/m)
        j_lo = j

        # Find next greatest Ys point k s.t. Ys[k] < X[i+1]
        k = _rank(Ys[j:], Xs[min(i+1, n)]) + j
        test_hi = np.abs(i/n - k/m)
        j_hi = k

        # Take the maximum distance, and corresponding index
        Tv[i] = np.max((test_lo, test_hi))
        js[i] = j_lo if np.argmax((test_lo, test_hi)) == 0 else j_hi

    return Tv, js


def _rank(A, k):
    """Return the number of keys in `A` strictly less than `k`."""
    assert all(A == sorted(A))
    lo = 0
    hi = len(A) - 1
    while lo <= hi:
        mid = (hi + lo) // 2
        if k < A[mid]:
            hi = mid - 1
        elif k > A[mid]:
            lo = mid + 1
        else:  # k == A[mid]
            return mid
    return lo
```

An example two-sample KS-test is shown in
Figure [2](#fig:ks_test){reference-type="ref" reference="fig:ks_test"}.

![The empirical cdfs of two independent random samples from
$\mathcal{N}\left( 0, 1 \right)$ and
$\mathcal{N}\left( 0, 2 \right)$. The test statistic $T_{n,m}$ is
shown by the double arrow.](/assets/images/ks2samp/ks_test.pdf){#fig:ks_test width="90%"}

### The Null Hypothesis

::: {.prop}
**Proposition 3**. *If $H_0$ is true, then
$$T_{n,m} = \sup_{0 \le x \le 1} \left| \frac{1}{n}\sum_{i=1}^{n} \mathbbm{1}\!\left\{U_i \le x\right\}
- \frac{1}{m}\sum_{j=1}^{m} \mathbbm{1}\!\left\{V_j \le x\right\} \right|.$$*
:::

::: {.proof}
*Proof.* By [\[eq:F_n\]](#eq:F_n){reference-type="eqref"
reference="eq:F_n"} and [\[eq:G_m\]](#eq:G_m){reference-type="eqref"
reference="eq:G_m"}, $$\label{eq:Tnm_supt}
    T_{n,m} = \sup_{t \in \mathbb{R}} \left| \frac{1}{n}\sum_{i=1}^{n} \mathbbm{1}\!\left\{X_i \le t\right\} - \frac{1}{m}\sum_{j=1}^{m} \mathbbm{1}\!\left\{Y_j \le t\right\} \right|.$$
To show the proposition is true, we make a change of variable. Let
$$x = F(t).$$ Then, $$t \in \mathbb{R}\implies x \in [0, 1].$$ Since $F$
and $G$ are continuous and monotonically increasing, $$\begin{aligned}
{3}
    X_i \le t &\iff F(X_i) \le F(t) \\
              &\iff U_i \le x &\quad&\text{(definition)}.
  \end{aligned}$$ Similarly, $$\begin{aligned}
{3}
    Y_i \le t &\iff G(Y_i) \le G(t) \\
              &\iff G(Y_i) \le F(t) &\quad&\text{(under $H_0$)} \\
              &\iff V_i \le x &\quad&\text{(definition)}.
  \end{aligned}$$ Substitution of these expressions
into [\[eq:Tnm_supt\]](#eq:Tnm_supt){reference-type="eqref"
reference="eq:Tnm_supt"} completes the proof. ◻
:::

### The Joint Distribution of the Samples

::: {#prop:Tnm .prop}
**Proposition 4**. *If $H_0$ is true, the joint distribution of
$U_1, \dots, U_n, V_1, \dots, V_m$ $(n+m)$ random variables is uniform
on $[0, 1]$.*
:::

::: {.proof}
*Proof.* $$\begin{aligned}
{3}
    \mathbb{P}\left[U_i \le t\right] &= \mathbb{P}\left[F(X_i) \le t\right] \\
                     &= \mathbb{P}\left[F(X_1) \le t\right] &\quad&\text{(i.i.d.)} \\
                     &= \mathbb{P}\left[G(X_1) \le t\right] &\quad&\text{(under $H_0$)} \\
                     &= \mathbb{P}\left[G(Y_1) \le t\right] &\quad&\text{(i.i.d.)} \\
                     &= \mathbb{P}\left[V_1 \le t\right] &\quad&\text{(definition)} \\
    \intertext{These probabilities can be rearranged to find the cdfs of $U$ and $V$}
                     &= \mathbb{P}\left[X_1 \le F^{-1}(t)\right] \\
                     &= F(F^{-1}(t)) &\quad&\text{(definition of cdf)} \\
                     &= t \\
    \therefore F_U(t) &= G_V(t) = t \\
    \implies f_{U,V}(t) &= \mathcal{U}\left(\left[ 0, 1 \right]\right) \tag*{\qedhere}
  \end{aligned}$$ ◻
:::

### The Test Statistic is Pivotal

Since Proposition [Proposition 4](#prop:Tnm){reference-type="ref"
reference="prop:Tnm"} has been shown to be true under the null
hypothesis $H_0$, and the distributions of $U_i$ and $V_j$ have been
shown to be $\mathcal{U}\left(\left[ 0, 1 \right]\right)$
independent of the distributions of the underlying samples $X_i$, $Y_j$,
we conclude that $T_{n,m}$ is *pivotal*, *i.e. *it does not itself
depend on the unknown distributions of the samples.

### Quantiles of the Test Statistic

Let $\alpha \in (0, 1)$ and $q_\alpha$ be the $(1 - \alpha)$-quantile of
the distribution of $T_{n,m}$ under $H_0$. The quantile $q_\alpha$ is
given by $$\begin{aligned}
  q_\alpha &= F^{-1}(1-\alpha) \\
           &= \inf\{x \colon F(x) \ge 1 - \alpha\} \\
           &\approx \min\{x \colon F_n(x) \ge 1 - \alpha\}, \quad n < \infty \\
           \implies q_\alpha \approx \hat{q}_\alpha &= \min_i \left\{
               T_{n,m}^{(i)} \colon \tfrac{i}{M} \ge 1 - \alpha \right\}\end{aligned}$$
where $M \in \mathbb{N}$ is large, and $T_{n,m}^{(i)}$ is the
$i^\text{th}$ value in a sorted sample of $M$ test statistics. Thus,
$q_\alpha$ can be approximated by choosing $i = \ceil{M(1 - \alpha)}$.
An algorithm to approximate $q_\alpha$ given $\alpha$ is as follows.

$n = \dim X$. $m = \dim Y$. $M \in \mathbb{N}$. $\alpha
    \in (0, 1)$. $q_\alpha \in [0, 1]$. $T_v \gets$ empty array of size
$n$ $X_s \gets$ sample of size $n$ from
$\mathcal{N}\left( 0, 1 \right)$. $Y_s \gets$ sample of size $m$ from
$\mathcal{N}\left( 0, 1 \right)$. $T_v^{(i)} \gets$ $T_{vs} \gets$
$j \gets \ceil*{M(1 - \alpha)}$ $T_{vs}^{(j)}$

A plot of the distribution of
$$\frac{T_{n,m}^M - \overline{T}_{n,m}^M}{\sqrt{\operatorname{Var}\left(T_{n,m}^M\right)}}$$
is shown in Figure [3](#fig:Tnm){reference-type="ref"
reference="fig:Tnm"} in comparison to a standard normal. The test
statistic distribution is skewed to the left, and has a longer right
tail than the standard normal. Since the asymptotic distribution of the
test statistic is not readily found in theory, we rely on simulation via
Algorithm [\[alg:ks_q\]](#alg:ks_q){reference-type="ref"
reference="alg:ks_q"} to estimate the quantiles.

![Empirical distribution of samples of the test statistic
$T_{n,m}$.](/assets/images/ks2samp/ks_dist.pdf){#fig:Tnm width="95%"}

### The Hypothesis Test

Given the aproximation for $\hat{q}_\alpha$ for $q_\alpha$ from
Algorithm [\[alg:ks_q\]](#alg:ks_q){reference-type="ref"
reference="alg:ks_q"}, we define a test with non-asymptotic level
$\alpha$ for $H_0$ vs. $H_1$:
$$\delta_\alpha = \mathbbm{1}\!\left\{T_{n,m} > \hat{q}_\alpha^{(n, M)}\right\}$$
where $T_{n,m}$ is found by
Algorithm [\[alg:ks_stat\]](#alg:ks_stat){reference-type="ref"
reference="alg:ks_stat"}. The p-value for this test is $$\begin{aligned}
  \text{p-value} &\coloneqq \mathbb{P}\left[Z \ge T_{n,m}\right] \\
  &\approx \frac{\#\{j = 1, \dots, M \colon T_{n,m}^{(j)} \ge T_{n,m}\}}{M}\end{aligned}$$
where $Z$ is a random variable distributed as $T_{n,m}$.

