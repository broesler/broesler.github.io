---
layout: post
title:  "Kolmogorov-Smirnov Test for Two Samples"
date: 2021-01-27
categories: statistics
tags: statistics hypothesis-testing python
---

<div style="visibility: hidden">
$$
\begin{align*}
\newcommand{\coloneqq}{\mathrel{\vcenter{:}}=}
\newcommand{\ceil}[1]{\left\lceil #1 \right\rceil}
\newcommand{\floor}[1]{\left\lfloor #1 \right\rfloor}
\end{align*}
$$
</div>


Consider two independent samples $X_1, \dots, X_n$, and
$Y_1, \dots, Y_m$ of independent, real-valued, continuous random
variables, and assume that the $X_i$’s are i.i.d. with some cdf $F$ and
that the $Y_i$’s are i.i.d. with some cdf $G$.[^1] We want to test
whether $F = G$. Consider the following hypotheses:

$$
\begin{align}
  H_0 \colon ``F = G" \nonumber\\
  H_1 \colon ``F \ne G"\nonumber
\end{align}
$$

 For simplicity, we will assume
that $F$ and $G$ are continuous and increasing.

## Example Experiment

An example experiment in which testing if two samples are from the same
distribution is of interest may be encountered in a lab setting where we
have two devices for measurement, and wish to determine if the errors
have the same distribution for our analysis.

## CDF Distributions

Let

$$
\begin{align}
  U_i &= F(X_i), \quad \forall i = 1, \dots, n, \nonumber\\
  V_j &= G(Y_j), \quad \forall j = 1, \dots, n.\nonumber
\end{align}
$$

<div class="prop" markdown=1>

**Proposition 1**. *The distribution of the cdf of a continuous random
variable is uniform on $[0,
  1]$.*

</div>

<div class="proof" markdown=1>

*Proof.* The distributions of $U_i$ and $V_j$ can be determined by
finding their cdfs. The cdf of $U_i$ is defined by
$F_U(t) \coloneqq \mathbb{P}\left[U_i \le t\right]$. Assuming that
$F(X)$ and $G(Y)$ are invertible, it follows that

$$
\begin{align}
  \mathbb{P}\left[U_i \le t\right] &= \mathbb{P}\left[F(X_i) \le t\right] &\quad&\text{(definition of $U_i$)} \nonumber\\
                   &= \mathbb{P}\left[X_i \le F^{-1}(t)\right] \nonumber\\
                   &= F(F^{-1}(t)) &\quad&\text{(definition of cdf)} \nonumber\\
                   &= t \nonumber\\
  \therefore F_U(t) &= t \nonumber\\
  \implies f_U(t) &= \mathcal{U}\left(\left[ 0, 1 \right]\right) \tag*{◻}
\end{align}
$$

</div>

Likewise, $f_V(t) = \mathcal{U}\left(\left[ 0, 1 \right]\right)$.

## Empirical CDFs

Let $F_n$ be the empirical cdf of $\{X_1, \dots, X_n\}$ and $G_m$ be the
empirical cdf of $\{Y_1, \dots, Y_m\}$.

### The Test Statistic

Let

$$
T_{n,m} = \sup_{t \in \mathbb{R}} \left| F_n(t) - G_m(t) \right|
$$

<div class="prop" markdown=1>

**Proposition 2**. *The test statistic $T_{n,m}$ can be written as the
maximum value of a finite set of numbers.*

</div>

<div class="proof" markdown=1>

*Proof.* By definition, the cdf

$$
\begin{align}
    F(t) &= \mathbb{P}\left[X \le t\right] \quad \forall t \in \mathbb{R}\nonumber\\
         &= \mathbb{E}\left[\mathbb{1}\!\left\{X \le t\right\}\right]. \nonumber\\
    \end{align}
$$

 By the Law of Large Numbers, the expectation can be approximated by the sample average, so we can define the *empirical cdf* as

$$
\begin{align}
    F_n(t) &= \frac{1}{n}\sum_{i=1}^{n} \mathbb{1}\!\left\{X_i \le t\right\}  \label{eq:F_n}
    \end{align}
$$

 Likewise,

$$
\begin{align}
    G_m(t) &= \frac{1}{m}\sum_{j=1}^{m} \mathbb{1}\!\left\{Y_j \le t\right\}.  \label{eq:G_m}
  \end{align}
$$

$$
\therefore T_{n,m} = \sup_{t \in \mathbb{R}} \left| \frac{1}{n}\sum_{i=1}^{n} \mathbb{1}\!\left\{X_i \le t\right\} - \frac{1}{m}\sum_{j=1}^{m} \mathbb{1}\!\left\{Y_j \le t\right\} \right|.
$$

The empirical cdfs $\eqref{eq:F_n}$ and $\eqref{eq:G_m}$ can also be
written

$$
\begin{align}
    F_n(t) &= \#\{i=1, \dots, n \colon X_i \le t\} \cdot \frac{1}{n} \nonumber\\
    G_m(t) &= \#\{i=1, \dots, m \colon Y_j \le t\} \cdot \frac{1}{m},\nonumber
  \end{align}
$$

 so the only values that the empirical cdfs can take
are the discrete sets

$$
\begin{align}
    F_n(i) &= \frac{i}{n} \quad \forall i = 1, \dots, n \nonumber\\
    G_m(j) &= \frac{j}{m} \quad \forall j = 1, \dots, m.\nonumber
  \end{align}
$$

 Therefore, the test statistic can be rewritten as the
maximum value of a finite set of numbers:

$$
\begin{split}
      T_{n,m} = \max_{i=0,\dots,n} \Bigg[
      &\max_{j=0,\dots,m} \left| \frac{i}{n} - \frac{j}{m} \right|
        \mathbb{1}\!\left\{Y^{(j)} \le X^{(i)} < Y^{(j+1)}\right\}, \\
      &\max_{k=j+1, \dots, m} \left| \frac{i}{n} - \frac{k}{m} \right|
        \mathbb{1}\!\left\{Y^{(k)} \le X^{(i+1)}\right\} \Bigg]
    \end{split}
$$

 where $X^{(i)}$ is the $i^\text{th}$ value in the
ordered set of data $X^{(1)} \le \cdots \le X^{(n)}$. The values
$X^{(0)}, Y^{(0)} \coloneqq -\infty$ are prepended to the otherwise
finite realizations to simplify the computation. <span class="qed_symbol">◻</span>

</div>

The following algorithm calculates the KS test statistic for two given
samples.

<div class="algorithm" id="alg:ks_stat">
<div class="alg_caption_div">
<span class="alg_title">Algorithm 1. </span>
<span class="alg_caption">Calculate the KS test statistic $T_{n,m}$ for two samples.</span>
</div>
<div class="algorithmic" markdown=1>
<p><span class="alg_command">Require:</span> $X, Y$ are vectors of real numbers.</p>
<p><span class="alg_command">Ensure:</span> $0 \le T_{n,m} \le 1$.</p>
<p><span class="alg_command">Procedure</span> <span class="alg_proc">KS2Sample</span>($X, Y$)</p>
<p>	$X_s \gets \{-\infty,$ <span class="alg_call">Sort</span>($X$)$\}$</p>
<p>	$Y_s \gets$ <span class="alg_call">Sort</span>($Y$)</p>
<p>	$n \gets \dim X_s$</p>
<p>	$m \gets \dim Y_s$</p>
<p>	$T_v \gets$ empty array of size $n$</p>
<p>	<span class="alg_command">for all</span> $i \in \{0, \dots, n\}$ <span class="alg_command">do</span></p>
<p>		$j \gets j$ + <span class="alg_call">Rank</span>($\{Y_s^{(\ell)}\}_{\ell=j}^m, X_s^{(i)}$) <span style="float: right">&#x25B7; Only search remaining $j$ values</span></p>
<p>		$k \gets j$ + <span class="alg_call">Rank</span>($\{Y_s^{(\ell)}\}_{\ell=j}^m, X_s^{(\min(i+1, n))}$)</p>
<p>		$\displaystyle{T_v^{(i)} \gets \max\left(\left|\frac{i}{n} - \frac{j}{m}\right|, \left|\frac{i}{n} - \frac{k}{m}\right|\right)}$</p>
<p>	<span class="alg_command">end for</span></p>
<p>	<span class="alg_command">Return</span> $\max_i T_v$</p>
<p><span class="alg_command">end procedure</span></p>
<p><span class="alg_command">Function</span> <span class="alg_proc">Rank</span>($A, k$)</p>
	
<p>	<span class="alg_command">Assert</span> $A$ is sorted in ascending order.</p>
<p>	<span class="alg_command">Return</span> $\#\{i=1,\dots,\dim A \colon k < A_i\}$</p>
<p><span class="alg_command">end function</span></p>
</div>
</div>

The following subroutine is an implementation of
<a href="#alg:ks_stat" data-reference-type="ref" data-reference="alg:ks_stat">Algorithm&nbsp;1</a>.
It computes an array of values $T_v(i)$ for each value of $X_i$. The
test statistic $T_{n,m}$ is the maximum of these values.

``` python
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
<a href="#fig:ks_test" data-reference-type="ref" data-reference="fig:ks_test">Figure 1</a>.

<figure>
<img src="/assets/images/ks2samp/ks_test.pdf" id="fig:ks_test" style="width:90.0%" /><figcaption><span class="fig_number">Figure 1</span>. The empirical cdfs of two independent random samples from <span class="math inline">\(\mathcal{N}\left( 0, 1 \right)\)</span> and <span class="math inline">\(\mathcal{N}\left( 0, 2 \right)\)</span>. The test statistic <span class="math inline">\(T_{n,m}\)</span> is shown by the double arrow.</figcaption>
</figure>

### The Null Hypothesis

<div class="prop" markdown=1>

**Proposition 3**. *If $H_0$ is true, then the test statistic

$$
T_{n,m} = \sup_{0 \le x \le 1} \left|
        \frac{1}{n}\sum_{i=1}^{n} \mathbb{1}\!\left\{U_i \le x\right\}
      - \frac{1}{m}\sum_{j=1}^{m} \mathbb{1}\!\left\{V_j \le x\right\}
    \right|.
$$

*

</div>

<div class="proof" markdown=1>

*Proof.* By $\eqref{eq:F_n}$ and $\eqref{eq:G_m}$,

$$
\begin{equation} \label{eq:Tnm_supt}
    T_{n,m} = \sup_{t \in \mathbb{R}} \left| \frac{1}{n}\sum_{i=1}^{n} \mathbb{1}\!\left\{X_i \le t\right\} - \frac{1}{m}\sum_{j=1}^{m} \mathbb{1}\!\left\{Y_j \le t\right\} \right|.
  \end{equation}
$$

To show the proposition is true, we make a change of variable. Let

$$
x = F(t).
$$

 Then,

$$
t \in \mathbb{R}\implies x \in [0, 1].
$$

 Since $F$
and $G$ are continuous and monotonically increasing,

$$
\begin{align}
    X_i \le t &\iff F(X_i) \le F(t) \nonumber\\
              &\iff U_i \le x &\quad&\text{(definition)}.\nonumber
  \end{align}
$$

 Similarly,

$$
\begin{align}
    Y_i \le t &\iff G(Y_i) \le G(t) \nonumber\\
              &\iff G(Y_i) \le F(t) &\quad&\text{(under $H_0$)} \nonumber\\
              &\iff V_i \le x &\quad&\text{(definition)}.\nonumber
  \end{align}
$$

 Substitution of these expressions
into $\eqref{eq:Tnm_supt}$ completes the proof. <span class="qed_symbol">◻</span>

</div>

### The Joint Distribution of the Samples

<div id="prop:Tnm" class="prop" markdown=1>

**Proposition 4**. *If $H_0$ is true, the joint distribution of
$U_1, \dots, U_n, V_1, \dots, V_m$ $(n+m)$ random variables is uniform
on $[0, 1]$.*

</div>

<div class="proof" markdown=1>

*Proof.*

$$
\begin{align}
    \mathbb{P}\left[U_i \le t\right] &= \mathbb{P}\left[F(X_i) \le t\right] \nonumber\\
                     &= \mathbb{P}\left[F(X_1) \le t\right] &\quad&\text{(i.i.d.)} \nonumber\\
                     &= \mathbb{P}\left[G(X_1) \le t\right] &\quad&\text{(under $H_0$)} \nonumber\\
                     &= \mathbb{P}\left[G(Y_1) \le t\right] &\quad&\text{(i.i.d.)} \nonumber\\
                     &= \mathbb{P}\left[V_1 \le t\right] &\quad&\text{(definition)} \nonumber\\
    \end{align}
$$

 These probabilities can be rearranged to find the cdfs of $U$ and $V$

$$
\begin{align}
                     &= \mathbb{P}\left[X_1 \le F^{-1}(t)\right] \nonumber\\
                     &= F(F^{-1}(t)) &\quad&\text{(definition of cdf)} \nonumber\\
                     &= t \nonumber\\
    \therefore F_U(t) &= G_V(t) = t \nonumber\\
    \implies f_{U,V}(t) &= \mathcal{U}\left(\left[ 0, 1 \right]\right) \tag*{◻}
  \end{align}
$$

</div>

### The Test Statistic is Pivotal

Since
<a href="#prop:Tnm" data-reference-type="ref" data-reference="prop:Tnm">Proposition 4</a>
has been shown to be true under the null hypothesis $H_0$, and the
distributions of $U_i$ and $V_j$ have been shown to be
$\mathcal{U}\left(\left[ 0, 1 \right]\right)$ independent of the
distributions of the underlying samples $X_i$, $Y_j$, we conclude that
$T_{n,m}$ is *pivotal*, *i.e. *it does not itself depend on the unknown
distributions of the samples.

### Quantiles of the Test Statistic

Let $\alpha \in (0, 1)$ and $q_\alpha$ be the $(1 - \alpha)$-quantile of
the distribution of $T_{n,m}$ under $H_0$. The quantile $q_\alpha$ is
given by

$$
\begin{align}
  q_\alpha &= F^{-1}(1-\alpha) \nonumber\\
           &= \inf\{x \colon F(x) \ge 1 - \alpha\} \nonumber\\
           &\approx \min\{x \colon F_n(x) \ge 1 - \alpha\}, \quad n < \infty \nonumber\\
           \implies q_\alpha \approx \hat{q}_\alpha &= \min_i \left\{\nonumber
               T_{n,m}^{(i)} \colon \tfrac{i}{M} \ge 1 - \alpha \right\}\nonumber
\end{align}
$$

where $M \in \mathbb{N}$ is large, and $T_{n,m}^{(i)}$ is the
$i^\text{th}$ value in a sorted sample of $M$ test statistics. Thus,
$q_\alpha$ can be approximated by choosing $i = \ceil{M(1 - \alpha)}$.
An algorithm to approximate $q_\alpha$ given $\alpha$ is as follows.

<div class="algorithm" id="alg:ks_q">
<div class="alg_caption_div">
<span class="alg_title">Algorithm 2. </span>
<span class="alg_caption">Approximate $q_\alpha$, the $(1 - \alpha)$-quantile of the distribution of $T_{n,m}$ under $H_0$.</span>
</div>
<div class="algorithmic" markdown=1>
<p><span class="alg_command">Require:</span> $n = \dim X$, $m = \dim Y$, $M \in \mathbb{N}$, and $\alpha \in (0, 1)$.</p>
<p><span class="alg_command">Ensure:</span> $q_\alpha \in [0, 1]$.</p>
<p><span class="alg_command">Procedure</span> <span class="alg_proc">KSQuantile</span>($n, m, M, \alpha$)</p>
<p>	$T_v \gets$ empty array of size $n$</p>
<p>	<span class="alg_command">for all</span> $i \in \{0,\dots,M\}$ <span class="alg_command">do</span></p>
<p>		$X_s \gets$ sample of size $n$ from $\mathcal{N}\left( 0, 1 \right)$.</p>
<p>		$Y_s \gets$ sample of size $m$ from $\mathcal{N}\left( 0, 1 \right)$.</p>
<p>		$T_v^{(i)} \gets$ <span class="alg_call">KS2Sample</span>($X_s, Y_s$) <span style="float: right">&#x25B7; defined in <a href="#alg:ks_stat">Algorithm&nbsp;1</a></span></p>
<p>	<span class="alg_command">end for</span></p>
<p>	$T_{vs} \gets$ <span class="alg_call">Sort</span>($T_v$)</p>
<p>	$j \gets \ceil{M(1 - \alpha)}$</p>
<p>	<span class="alg_command">Return</span> $T_{vs}^{(j)}$</p>
<p><span class="alg_command">end procedure</span></p>
</div>
</div>

A plot of the distribution of

$$
\frac{T_{n,m}^M - \overline{T}_{n,m}^M}{\sqrt{\operatorname{Var}\left(T_{n,m}^M\right)}}
$$

is shown in
<a href="#fig:Tnm" data-reference-type="ref" data-reference="fig:Tnm">Figure 2</a>
in comparison to a standard normal. The test statistic distribution is
skewed to the left, and has a longer right tail than the standard
normal. Since the asymptotic distribution of the test statistic is not
readily found in theory, we rely on simulation via
<a href="#alg:ks_q" data-reference-type="ref" data-reference="alg:ks_q">Algorithm&nbsp;2</a>
to estimate the quantiles.

<figure>
<img src="/assets/images/ks2samp/ks_dist.pdf" id="fig:Tnm" style="width:95.0%" /><figcaption><span class="fig_number">Figure 2</span>. Empirical distribution of samples of the test statistic <span class="math inline">\(T_{n,m}\)</span>.</figcaption>
</figure>

### The Hypothesis Test

Given the aproximation for $\hat{q}_\alpha$ for $q_\alpha$ from
<a href="#alg:ks_q" data-reference-type="ref" data-reference="alg:ks_q">Algorithm&nbsp;2</a>,
we define a test with non-asymptotic level $\alpha$ for $H_0$ vs. $H_1$:

$$
\delta_\alpha = \mathbb{1}\!\left\{T_{n,m} > \hat{q}_\alpha^{(n, M)}\right\}
$$

where $T_{n,m}$ is found by
<a href="#alg:ks_stat" data-reference-type="ref" data-reference="alg:ks_stat">Algorithm&nbsp;1</a>.
The p-value for this test is

$$
\begin{align}
  \text{p-value} &\coloneqq \mathbb{P}\left[Z \ge T_{n,m}\right] \nonumber\\
  &\approx \frac{\#\{j = 1, \dots, M \colon T_{n,m}^{(j)} \ge T_{n,m}\}}{M}\nonumber
\end{align}
$$

where $Z$ is a random variable distributed as $T_{n,m}$.

[^1]: Note that the two samples may have different sizes (if $n \ne m$).
