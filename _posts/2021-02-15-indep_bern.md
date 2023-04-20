---
layout: post
title:  "Tests of Independence: Part 1 – Bernoulli Random Variables"
date:   "2021-02-15 16:21:56 -0500"
categories: statistics
tags: statistics hypothesis-testing python
reading_time: 29
---

<p class="message">
This post is the first in a series on tests for independence of random
variables.
</p> 

Let's start with the simplest case of Bernoulli random variables. I'll
explore continuous random variables in a later post. This post roughly follows
[MIT's 18.650][ocw], Problem Set 6, although I've written it in a 
proposition-proof style instead of more verbose prose.

## Problem Statement

Let $X, Y$ be two Bernoulli random variables, not necessarily
independent, and let

$$
\begin{align*}
  p &= \mathbb{P}\left[X=1\right], \\
  q &= \mathbb{P}\left[Y=1\right], \text{ and} \\
  r &= \mathbb{P}\left[X=1, Y=1\right].
\end{align*}
$$

We now look to
define a test to show that $X$ and $Y$ are independent.

## Condition for Independence

<div class="prop" markdown=1>

**Proposition 1**. *$X \indep Y \iff r = pq$.*

</div>

<div class="proof" markdown=1>

*Proof.* Two random variables are independent iff

$$
\begin{equation*}
    \mathbb{P}\left[X \cap Y\right] = \mathbb{P}\left[X\right]\mathbb{P}\left[Y\right] \label{eq:indep}
    % \iff \Prob{X} = \frac{\Prob{X,Y}}{\Prob{Y}} = \Prob{X | Y}
  \end{equation*}
$$

By the given definition of the Bernoulli random variables,

$$
\begin{align*}
    \mathbb{P}\left[X \cap Y\right] &\equiv \mathbb{P}\left[X, Y\right] \equiv \mathbb{P}\left[X=1, Y=1\right] = r \\
    \mathbb{P}\left[X\right] &\equiv \mathbb{P}\left[X=1\right] = p \\
    \mathbb{P}\left[Y\right] &\equiv \mathbb{P}\left[Y=1\right] = q \\
    \therefore r = pq &\iff \mathbb{P}\left[X,Y\right] = \mathbb{P}\left[X\right]\mathbb{P}\left[Y\right] \\
    &\iff X \indep Y  \tag*{◻}
  \end{align*}
$$

</div>

## Test for Independence

Let $(X_1, Y_1), \dots, (X_n, Y_n)$ be a sample of $n$ i.i.d. copies of
$(X, Y)$ (*i.e. *$X_i \indep X_j$ for $i \ne j$, but $X_i$ may not be
independent of $Y_i$). Based on this sample, we want to test whether
$X \indep Y$, *i.e. *whether $r = pq$.

### Estimators of the Parameters

Define the estimators:

$$
\begin{align*}
  \hat{p}&= \frac{1}{n}\sum_{i=1}^{n} X_i, \\
  \hat{q}&= \frac{1}{n}\sum_{i=1}^{n} Y_i, \\
  \hat{r}&= \frac{1}{n}\sum_{i=1}^{n} X_i Y_i.
\end{align*}
$$

<div class="prop" markdown=1>

**Proposition 2**. *These estimators $\hat{p}$, $\hat{q}$, and $\hat{r}$
are consistent estimators of the true values $p$, $q$, and $r$.*

</div>

<div class="proof" markdown=1>

*Proof.* To show that an estimator is *consistent*, we must prove that
it converges to the true value of the parameter in the limit as
$n \to \infty$. Since the sequence of $X_i$’s are i.i.d., we can use the
Weak Law of Large Numbers (LLN) to prove that $\hat{p}$ converges to
$p$.

<div id="eq:LLN" class="theorem" markdown=1>

**Theorem 1** (Weak Law of Large Numbers). *If the sequence of random
variables $X_1, \dots, X_n$ are i.i.d., then*

$$
\frac{1}{n}\sum_{i=1}^{n} X_i \xrightarrow[n \to \infty]{\mathbb{P}}\mathbb{E}\left[X\right].
$$

</div>

The expectation of $X$ is given by

$$
\begin{align*}
    \mathbb{E}\left[X\right] &= \mathbb{E}\left[\mathop{\mathrm{Ber}}(p)\right] &\quad&\text{(given)} \\
          &= p &\quad&\text{(definition of Bernoulli r.v.)} \\
    \therefore \frac{1}{n}\sum_{i=1}^{n} X_i &\xrightarrow[n \to \infty]{\mathbb{P}}p &\quad&\text{(LLN)} \\
    \implies \hat{p}&\xrightarrow[n \to \infty]{\mathbb{P}}p.
  \end{align*}
$$

Likewise
$\hat{q}\xrightarrow[n \to \infty]{\mathbb{P}}q$.

To show that $\hat{r}$ converges to $r$, let $R \coloneqq X Y$ be a
Bernoulli random variable with parameter
$r = \mathbb{P}\left[X=1, Y=1\right]$, so that the estimator

$$
\begin{equation*} \label{eq:rhat}
    \hat{r}= \frac{1}{n}\sum_{i=1}^{n} X_i Y_i = \frac{1}{n}\sum_{i=1}^{n} R_i.
  \end{equation*}
$$

Note that the values of $R_i$ *are* i.i.d. since each pair $(X_i, Y_i)$
are i.i.d., even though $X_i$ and $Y_j$ may not be independent for
$i \ne j$. As before, we apply the Law of Large Numbers to the average
of $R_i$’s. The expectation of $R$ is

$$
\begin{align*}
    \mathbb{E}\left[R\right] &= \mathbb{E}\left[\mathop{\mathrm{Ber}}(r)\right] &\quad&\text{(definition)} \\
          &= r &\quad&\text{(definition of Bernoulli r.v.)} \\
    \therefore \frac{1}{n}\sum_{i=1}^{n} R_i &\xrightarrow[n \to \infty]{\mathbb{P}}r &\quad&\text{(LLN)} \\
    \implies \hat{r}&\xrightarrow[n \to \infty]{\mathbb{P}}r.
  \end{align*}
$$

Thus, each estimator $(\hat{p}, \hat{q}, \hat{r})$
converges to its respective parameter $(p, q, r)$. <span class="qed_symbol">◻</span>

</div>

### Asymptotic Normality of the Estimators

<div id="prop:normality_pqr" class="prop" markdown=1>

**Proposition 3**. *The vector of estimators
$(\hat{p}, \hat{q}, \hat{r})$ is asymptotically normal,* i.e.

$$
\sqrt{n} ((\hat{p}, \hat{q}, \hat{r}) - (p, q, r))
        \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, \operatorname{Cov}\left((\hat{p}, \hat{q}, \hat{r})\right) \right).
$$

</div>

*Proof.* To prove that the vector of estimators is asymptotically
normal, we employ the Central Limit Theorem (CLT).

<div class="theorem" markdown=1>

**Theorem 2** (Central Limit Theorem). *Let $X_1, \dots, X_n$ be a
sequence of i.i.d. random vectors $X_i \in \mathbb{R}^k$, and
$\overline{X}_n= \frac{1}{n}\sum_{i=1}^{n} X_i$. Then,*

$$
\sqrt{n}(\overline{X}_n- \mathbb{E}\left[X\right]) \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, \Sigma \right).
$$

*where $\Sigma$ is the $k$-by-$k$ matrix
$\operatorname{Cov}\left(X\right)$.*

</div>

By the CLT,

$$
\begin{equation*} \label{eq:CLT}
  \sqrt{n} ((\hat{p}, \hat{q}, \hat{r}) - \mathbb{E}\left[(\hat{p}, \hat{q}, \hat{r})\right]) \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, \Sigma \right)
\end{equation*}
$$

where $\Sigma$ is the 3-by-3 symmetric covariance matrix, defined as

$$
\begin{equation*}
  \Sigma \coloneqq
  \begin{bmatrix}
    \operatorname{Var}\left(\hat{p}\right) & \operatorname{Cov}\left(\hat{p}, \hat{q}\right) & \operatorname{Cov}\left(\hat{p}, \hat{r}\right) \\
    \cdot & \operatorname{Var}\left(\hat{q}\right) & \operatorname{Cov}\left(\hat{q}, \hat{r}\right) \\
    \cdot & \cdot & \operatorname{Var}\left(\hat{r}\right)
  \end{bmatrix}.
\end{equation*}
$$

We first need to determine the expectations of the estimators.

<div class="prop" markdown=1>

**Proposition 4**. *The expectation of the estimator $\hat{p}$ is
$\mathbb{E}\left[\hat{p}\right] = p$.*

</div>

<div class="proof" markdown=1>

*Proof.*

$$
\begin{align*}
    \mathbb{E}\left[\hat{p}\right] &= \mathbb{E}\left[\frac{1}{n}\sum_{i=1}^{n} X_i\right] &\quad&\text{(definition)} \\
              &= \frac{1}{n} \sum_{i=1}^{n} \mathbb{E}\left[X_i\right] &\quad&\text{(linearity of expectation)} \\
              &= \frac{1}{n} n \mathbb{E}\left[X\right] &\quad&\text{(i.i.d.)} \\
              &= \mathbb{E}\left[\mathop{\mathrm{Ber}}(p)\right] &\quad&\text{(definition)} \\
    \implies \mathbb{E}\left[\hat{p}\right] &= p &\quad&\text{(definition)} \tag*{◻}
  \end{align*}
$$

</div>

Similarly, $\mathbb{E}\left[\hat{q}\right] = q$ and
$\mathbb{E}\left[\hat{r}\right] = r$. This proposition also shows that
the estimators are *unbiased*, since
$\mathbb{E}\left[\hat{p}- p\right] = 0$, *etc*.

We now determine the entries in the covariance matrix to complete the
proof of asymptotic normality.

<div class="prop" markdown=1>

**Proposition 5**. *The variance of $\hat{p}$ is given by
$\operatorname{Var}\left(\hat{p}\right) = \frac{1}{n} p(1 - p)$.*

</div>

<div class="proof" markdown=1>

*Proof.* Using the definition of $\hat{p}$,

$$
\begin{align*}
    \operatorname{Var}\left(\hat{p}\right) &= \operatorname{Var}\left(\frac{1}{n}\sum_{i=1}^{n} X_i\right) &\quad&\text{(definition)} \\
                &= \frac{1}{n^2}\operatorname{Var}\left(\sum_{i=1}^{n} X_i\right) &\quad&\text{(variance rule)} \\
                &= \frac{1}{n^2}\sum_{i=1}^{n}\operatorname{Var}\left(X_i\right) &\quad&\text{(i.i.d.)} \\
                &= \frac{1}{n^2}n\operatorname{Var}\left(X\right) &\quad&\text{(i.i.d.)} \\
   \therefore \operatorname{Var}\left(\hat{p}\right) &= \frac{1}{n} p(1-p) &\quad&\text{(variance of $\mathop{\mathrm{Ber}}(p)$)}
  \end{align*}
$$

Likewise,
$\operatorname{Var}\left(\hat{q}\right) = \frac{1}{n} q(1 - q)$, and
$\operatorname{Var}\left(\hat{r}\right) = \frac{1}{n} r(1 - r)$. <span class="qed_symbol">◻</span>

</div>

<div class="prop" markdown=1>

**Proposition 6**. *The covariance of $\hat{p}$ and $\hat{q}$ is given
by $\operatorname{Cov}\left(\hat{p}, \hat{q}\right) = r - pq$.*

</div>

<div class="proof" markdown=1>

*Proof.*

$$
\begin{align*}
    \operatorname{Cov}\left(\hat{p}, \hat{q}\right) &= \operatorname{Cov}\left(\frac{1}{n}\sum_{i=1}^{n} X_i, \frac{1}{n}\sum_{i=1}^{n} Y_i\right) \\
    &= \frac{1}{n^2} \operatorname{Cov}\left(\sum_{i=1}^{n} X_i, \sum_{i=1}^{n} Y_i\right) &\quad&\text{(covariance property)} \\
    &= \frac{1}{n^2} \sum_{i=1}^{n} \sum_{j=1}^{n} \operatorname{Cov}\left(X_i, Y_j\right) &\quad&\text{(bilinearity of covariance)} \\
    &= \frac{1}{n^2} n^2 \operatorname{Cov}\left(X, Y\right) &\quad&\text{(identically distributed)} \\
    &= \operatorname{Cov}\left(X, Y\right)  \\
    &= \mathbb{E}\left[XY\right] - \mathbb{E}\left[X\right]\mathbb{E}\left[Y\right] &\quad&\text{(definition of covariance)} \\
    &= \mathbb{E}\left[R\right] - \mathbb{E}\left[X\right]\mathbb{E}\left[Y\right] &\quad&\text{(definition of $R$)} \\
    \therefore \operatorname{Cov}\left(\hat{p}, \hat{q}\right) &= r - pq \tag*{◻}
  \end{align*}
$$

</div>

<div class="prop" markdown=1>

**Proposition 7**. *The covariance of $\hat{p}$ and $\hat{r}$ is given
by $\operatorname{Cov}\left(\hat{p}, \hat{r}\right) = r(1 - p)$.*

</div>

<div class="proof" markdown=1>

*Proof.*

$$
\begin{align*}
    \operatorname{Cov}\left(\hat{p}, \hat{r}\right) &= \operatorname{Cov}\left(\frac{1}{n}\sum_{i=1}^{n} X_i, \frac{1}{n}\sum_{i=1}^{n} R_i\right) \\
    &= \frac{1}{n^2} \operatorname{Cov}\left(\sum_{i=1}^{n} X_i, \sum_{i=1}^{n} R_i\right) &\quad&\text{(covariance property)} \\
    &= \frac{1}{n^2} \sum_{i=1}^{n} \sum_{j=1}^{n} \operatorname{Cov}\left(X_i, R_j\right) &\quad&\text{(bilinearity of covariance)} \\
    &= \frac{1}{n^2} n^2 \operatorname{Cov}\left(X, R\right) &\quad&\text{(identically distributed)} \\
    &= \operatorname{Cov}\left(X, R\right)  \\
    &= \mathbb{E}\left[X R\right] - \mathbb{E}\left[X\right]\mathbb{E}\left[R\right] &\quad&\text{(definition of covariance)} \\
    &= \mathbb{E}\left[X R\right] - pr &\quad&\text{(given)} \\
    &= \mathbb{E}\left[X (X Y)\right] - pr &\quad&\text{(definition of $R$)} \\
    \end{align*}
$$

Since $X \sim \mathop{\mathrm{Ber}}(p) \in \\\{0, 1\\\}$, $X^2 = X$, so we have

$$
\begin{align*}
    &= \mathbb{E}\left[X Y\right] - pr \\
    &= r - pr \\
    \therefore \operatorname{Cov}\left(\hat{p}, \hat{r}\right) &= r(1 - p) \tag*{◻}
  \end{align*}
$$

</div>

Similarly, $\operatorname{Cov}\left(\hat{q}, \hat{r}\right) = r(1 - q)$.

The entire asymptotic covariance matrix is then

$$
\begin{equation} \label{1}\tag{1}
  \Sigma =
  \begin{bmatrix}
    p(1-p) & r - pq & r(1-p) \\
    \cdot & q(1-q) & r(1-q) \\
    \cdot & \cdot & r(1-r)
  \end{bmatrix}.
\end{equation}
$$

Since we have determined the expectation
$\mathbb{E}\left[(\hat{p}, \hat{q}, \hat{r})\right] = (p, q,
r)$, and the covariance matrix $\Sigma$ in terms of $p$, $q$, and $r$,
we conclude that
<a href="#prop:normality_pqr" data-reference-type="ref" data-reference="prop:normality_pqr">Proposition 3</a>
is true, and the vector of estimators $(\hat{p}, \hat{q}, \hat{r})$ is
asymptotically normal. <span class="qed_symbol">■</span>

### The Delta Method

<div id="prop:delta" class="prop" markdown=1>

**Proposition 8**.
$$
\sqrt{n}\left((\hat{r}- \hat{p}\hat{q}) - (r - pq)\right) \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, V \right)
$$

*where $V$ depends only on $p$, $q$, and $r$.*

</div>

<div class="proof" markdown=1>

*Proof.* Let $\hat{\theta}$ and $\theta$ be vectors in $\mathbb{R}^3$

$$
\hat{\theta} = \begin{bmatrix} \hat{p}\\ \hat{q}\\ \hat{r}
\end{bmatrix} \text{, and }
    \theta = \begin{bmatrix} p \\ q \\ r \end{bmatrix}.
$$

From our proof
of
<a href="#prop:normality_pqr" data-reference-type="ref" data-reference="prop:normality_pqr">Proposition 3</a>,
we have

$$
\begin{align*}
    \sqrt{n}(\hat{\theta} - \theta) &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, \Sigma \right) &\quad&\text{(CLT)} \\
    \implies \sqrt{n}(g(\hat{\theta}) - g(\theta)) &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, \nabla g(\theta)^\top\Sigma
      \nabla g(\theta) \right) &\quad&\text{(Delta method)}
  \end{align*}
$$

for any differentiable function
$g \colon \mathbb{R}^k \to \mathbb{R}$, and $\Sigma$ given by
Equation $\eqref{1}$. Define the function

$$
\begin{equation} \label{2}\tag{2}
    g(u, v, w) = w - uv
  \end{equation}
$$

such that

$$
\begin{align*}
    g(\hat{\theta}) &= \hat{r}- \hat{p}\hat{q}, \\
    g(\theta) &= r - pq.
  \end{align*}
$$

The gradient of $g(\theta)$ is then

$$
\nabla g(u,v,w) = \begin{bmatrix} -v \\ -u \\ 1 \end{bmatrix}
    \implies \nabla g(\theta) = \begin{bmatrix} -q \\ -p \\ 1 \end{bmatrix}
$$

The asymptotic variance
$V = \nabla g(\theta)^\top\Sigma \nabla g(\theta)$, which we will now
show is a function only of the parameters $(p, q, r)$.

$$
\begin{align*}
    V &= \begin{bmatrix} -q & -p & 1 \end{bmatrix}
    \begin{bmatrix}
      p(1-p) & r - pq & r(1-p) \\
      \cdot & q(1-q) & r(1-q) \\
      \cdot & \cdot & r(1-r)
    \end{bmatrix}
    \begin{bmatrix} -q \\ -p \\ 1 \end{bmatrix}  \label{eq:V_mat} \\
    &= \begin{bmatrix} -q & -p & 1 \end{bmatrix}
    \begin{bmatrix}
      -qp(1-p) - p(r - pq) + r(1-p) \\
      -q(r - pq) - pq(1-q) + r(1-q) \\
      -qr(1-p) - pr(1-q) + r(1-r)
    \end{bmatrix} \\
    &= \begin{bmatrix} -q & -p & 1 \end{bmatrix}
    \begin{bmatrix}
      (r - pq)(1 - 2p) \\
      (r - pq)(1 - 2q) \\
      r((1-p)(1-q) - (r-pq))
    \end{bmatrix} \\
    &= -q(r - pq)(1 - 2p) - p(r - pq)(1 - 2q))  \\
    &\,\quad + r((1-p)(1-q) - (r-pq))  \\
    \therefore V &= (r - pq)[-q(1 - 2p) - p(1 - 2q) - r] + r(1-p)(1-q) \label{3}\tag{3}
  \end{align*}
$$

which is a function only of $(p, q, r)$. <span class="qed_symbol">◻</span>

</div>

### The Null Hypothesis

Consider the hypotheses

$$
\begin{align*}
  H_0 \colon X \indep Y \\
  H_1 \colon X \nindep Y
\end{align*}
$$

<div id="prop:V_H0" class="prop" markdown=1>

**Proposition 9**. *If $H_0$ is true, then $V = pq(1-p)(1-q)$.*

</div>

<div class="proof" markdown=1>

*Proof.* Under $H_0$, $r = pq$. Using the previous expression for $V$,
Equation $\eqref{3}$, replace $r$ by $pq$ to find

$$
V = (pq - pq)[-q(1 - 2p) - p(1 - 2q) - pq] + pq(1-p)(1-q).
$$

The first
term is identically 0, so

$$
V = pq(1-p)(1-q). \tag*{◻}
$$

</div>

<div id="prop:Vhat" class="prop" markdown=1>

**Proposition 10**. *Given*

$$
V = pq(1-p)(1-q),
$$

*a consistent estimator
is given by*

$$
\hat{V} = \hat{p}\hat{q}(1 - \hat{p})(1 - \hat{q}).
$$

</div>

<div class="proof" markdown=1>

*Proof.* To prove that $\hat{V}$ converges to $V$, we employ the
Continuous Mapping Theorem.

<div class="theorem" markdown=1>

**Theorem 3** (Continuous Mapping Theorem). *Let $X \in \mathbb{R}^n$ be
a vector of random variables, and $g \colon \mathbb{R}^n \to \mathbb{R}$
be a continuous function. Let $X_n = X_1, X_2, \dots$ be a sequence of
random vectors. If $X_n \xrightarrow[n \to \infty]{\mathbb{P}}X$, then
$g(X_n) \xrightarrow[n \to \infty]{\mathbb{P}}g(X)$.*

</div>

Since $\hat{p}\xrightarrow[n \to \infty]{\mathbb{P}}p$ and
$\hat{q}\xrightarrow[n \to \infty]{\mathbb{P}}q$,
$\hat{V}(\hat{p}, \hat{q}) \xrightarrow[n \to \infty]{\mathbb{P}}V(p, q)$. <span class="qed_symbol">◻</span>

</div>

### A Hypothesis Test

<div class="prop" markdown=1>

**Proposition 11**. *Given $\alpha \in (0, 1)$, we propose the test
statistic*

$$
T_n \coloneqq \frac{\sqrt{n}(\hat{r}- \hat{p}\hat{q})}{\sqrt{\hat{V}}} \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right)
$$

*where $\hat{V}$ is given by
<a href="#prop:Vhat" data-reference-type="ref" data-reference="prop:Vhat">Proposition 10</a>,
and $t_{n-1}$ is Student’s $t$-distribution with $n-1$ degrees of
freedom.*

</div>

<div class="proof" markdown=1>

*Proof.*
<a href="#prop:delta" data-reference-type="ref" data-reference="prop:delta">Proposition 8</a>
gives the distribution of $g(\theta)$ (given by Equation $\eqref{2}$)
under $H_0$. Assume that $p, q \in (0, 1)$ s.t. $V > 0$.

$$
\begin{align}
  \sqrt{n}\left((\hat{r}- \hat{p}\hat{q}) - (r - pq)\right) &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, V \right) \nonumber \\
  \sqrt{n}(\hat{r}- \hat{p}\hat{q}) &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, V \right) &\quad&\text{($r = pq$ under $H_0$)} \nonumber \\
  \sqrt{n}\frac{(\hat{r}- \hat{p}\hat{q})}{\sqrt{V}} &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right)  \label{4}\tag{4}
\end{align}
$$

The asymptotic variance $V$, however, is unknown, so we divide the
estimator by $\sqrt{\frac{\hat{V}}{V}}$ to get an expression that will
evaluate to our desired test statistic

$$
T_n = \frac{\displaystyle \sqrt{n}\frac{(\hat{r}- \hat{p}\hat{q})}{\sqrt{V}}}{\displaystyle \sqrt{\frac{\hat{V}}{V}}}
$$

Given this expression, we can determine the distribution of $T_n$.
Equation $\eqref{4}$ shows that the numerator is a standard
normal random variable. [Cochran’s
theorem](https://en.wikipedia.org/wiki/Cochran's_theorem) gives the
distribution of the denominator.

<div class="lemma" markdown=1>

**Lemma 4** (Result of Cochran’s Theorem). *If $X_1, \dots, X_n$ are
i.i.d. random variables drawn from the distribution
$\mathcal{N}\left( \mu, \sigma^2 \right)$, and
$S_n^2 \coloneqq \sum_{i=1}^{n} (X_i - \overline{X}_n)^2$, then*

$$
\overline{X}_n\indep S_n,
$$

*and*

$$
\frac{n S_n^2}{\sigma^2} \sim \chi^2_{n-1}.
$$

</div>

Since $\hat{V}$ and $V$ describe the sample variance and variance of a
(asymptotically) normal distribution, $T_n$ is asymptotically
characterized by

$$
T_n \xrightarrow[n \to \infty]{(d)}\frac{\displaystyle \mathcal{N}\left( 0, 1 \right)}{\displaystyle \sqrt{\frac{\chi^2_{n-1}}{n}}}
$$

which is the definition of a random variable drawn from *Student’s
T-distribution* with $n-1$ degrees of freedom. In this case, however,
the normality of the underlying random variables is asymptotic, so the
$t_{n-1}$ distribution approaches a standard normal distribution

$$
\begin{align*}
  T_n &\xrightarrow[n \to \infty]{(d)}t_{n-1} \\
  t_{n-1} &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right)  \label{5}\tag{5} \\
  \implies T_n &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right) \tag*{◻}
\end{align*}
$$

A proof of Equation $\eqref{5}$ is given in the 
<a href="#app:A">Appendix</a>.

</div>

Given the test statistic $T_n$, define the rejection region

$$
R_\psi = \left\{ \hat{\theta} \colon |T_n| > q_{\alpha/2} \right\}
$$

where

$$
q_{\alpha/2} = \Phi^{-1}\left(1 - \frac{\alpha}{2}\right)
$$

is
the $\left(1-\frac{\alpha}{2}\right)$-quantile of the standard normal
$\mathcal{N}\left( 0, 1 \right)$ distribution.

We would like to know whether the facts of being happy and being in a
relationship are independent of each other. In a given population, 1000
people (aged at least 21 years old) are sampled and asked two questions:
“Do you consider yourself as happy?” and “Are you involved in a
relationship?”. The answers are summarized in
<a href="#tab:tab1" data-reference-type="ref" data-reference="tab:tab1">Table 1</a>.

<div id="tab:tab1" style="width: auto; margin: auto;" markdown=1>
<span class="marginnote-wrapper"><label for="sn-0" class="margin-toggle">&#8862;</label>
<input type="checkbox" id="sn-0" class="margin-toggle"/>
<span class="marginnote">Data sourced from [MIT's 18.650][ocw].</span></span>


|                           | **Happy** | **Not Happy** | **Total** |
|--------------------------:|:---------:|:-------------:|:---------:|
|     **In a Relationship** |    205    |      301      |    506    |
| **Not in a Relationship** |    179    |      315      |    494    |
|                 **Total** |    384    |      616      |   1000    |

</div>

The values of our estimators are as follows:

$$
\begin{align*}
  \hat{p}&= \frac{\text{\# Happy}}{N} = \frac{384}{1000} = 0.384 \\
  \hat{q}&= \frac{\text{\# In a Relationship}}{N} = \frac{506}{1000} = 0.506 \\
  \hat{r}&= \frac{\text{\# Happy} \cap \text{\# In a Relationship}}{N} = \frac{205}{1000} = 0.205.
\end{align*}
$$

The estimate of the asymptotic variance of the test statistic is

$$
\hat{V} = \hat{p}\hat{q}(1-\hat{p})(1-\hat{q}) = (0.384)(0.506)(1 - 0.384)(1 - 0.506)
  = 0.05913,
$$

giving the test statistic

$$
T_n = \frac{\sqrt{n}(\hat{r}- \hat{p}\hat{q})}{\sqrt{\hat{V}}}
  = \frac{\sqrt{1000}(0.205 - 0.384\cdot0.506)}{\sqrt{0.05913}} = 1.391.
$$

The standard normal quantile at 
$\alpha = 0.05$ is 
$q_{\alpha/2} = \Phi^{-1}\left(1 - \frac{\alpha}{2}\right) = 1.96$,
so the test result is

$$
|T_n| = 1.391 < q_{\alpha/2} = 1.96
$$

so we *fail to reject $H_0$ at
the 5% level*. The $p$-value of the test is

$$
\begin{align*}
  \text{$p$-value} &\coloneqq \mathbb{P}\left[Z > |T_n|\right] &\quad&\text{($Z \sim \mathcal{N}\left( 0, 1 \right)$)}  \\
                   &= \mathbb{P}\left[Z \le |T_n|\right] &\quad&\text{(symmetry)} \\
                   &= \mathbb{P}\left[Z \le -T_n\right] + \mathbb{P}\left[Z > T_n\right] \\
                   &= 2\Phi(-|T_n|) \\
  \implies \text{$p$-value} &= 0.1642.
\end{align*}
$$

In other words,
the lowest level at which we could reject the null hypothesis is at
$\alpha = \text{$p$-value} = 0.1642 = 16.42\%$.


<h2 id="app:A">Appendix A: Additional Proofs</h2>

<div class="app-wrapper">

<label for="app-A" class="app-toggle">Show/Hide Appendix &#8862;</label>
<input type="checkbox" id="app-A" class="app-toggle"/>


<div class="appendix" markdown=1>

<div class="prop" markdown=1>

**Proposition 12**. *A $t$-distribution with $n$ degrees of freedom
approaches a standard normal distribution as $n$ approaches infinity:*

$$
t_n \xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right).
$$

</div>

*Proof.* Student’s $t$-distribution with $\nu$ degrees of freedom is
defined as the distribution of the random variable $T$ such that

$$
t_{\nu} \sim T = \frac{\displaystyle Z}{\displaystyle \sqrt{\frac{V}{\nu}}}
$$

where $Z \sim \mathcal{N}\left( 0, 1 \right)$, $V \sim \chi^2_{\nu}$,
and $Z \indep V$.

Let $X_1, \dots, X_n \sim \mathcal{N}\left( \mu, \sigma^2 \right)$ be a
sequence of i.i.d. random variables. Define the sample mean and sample
variance

$$
\begin{align*}
  \overline{X}_n&\coloneqq \frac{1}{n}\sum_{i=1}^{n} X_i \\
  S_n^2 &\coloneqq \frac{1}{n}\sum_{i=1}^{n} (X_i - \overline{X}_n)^2.
\end{align*}
$$

Let the random variables

$$
\begin{align*}
  Z &= \frac{\sqrt{n}(\overline{X}_n- \mu)}{\sigma} \\
  V &= \frac{n S_n^2}{\sigma^2}.
\end{align*}
$$

such that
$Z \sim \mathcal{N}\left( 0, 1 \right)$ by the Central Limit Theorem,
and $V \sim
\chi^2_{n-1}$ by Cochran’s Theorem (which also shows that $Z \indep V$).
Then, the $t$-distribution is *pivotal*

$$
t_{n-1} = \frac{\displaystyle Z}{\displaystyle \sqrt{\frac{V}{n-1}}}.
$$

<div class="lemma" markdown=1>
**Lemma 5**. *The sample variance converges in probability to the variance,*

$$
S_n^2 \xrightarrow[n \to \infty]{\mathbb{P}}\sigma^2.
$$

</div>  

<div class="proof" markdown=1>
*Proof.*

$$
\begin{align*} 
S_n^2 &\coloneqq \frac{1}{n}\sum_{i=1}^{n}(X_i - \overline{X}_n)^2 \\
 &= \frac{1}{n}\sum_{i=1}^{n}(X_i^2 - 2 \overline{X}_nX_i + \overline{X}_n^2) \\
 &= \frac{1}{n}\sum_{i=1}^{n} X_i^2 - \frac{1}{n}\sum_{i=1}^{n} 2 \overline{X}_nX_i + \frac{1}{n}\sum_{i=1}^{n} \overline{X}_n^2  \\
 &= \frac{1}{n}\sum_{i=1}^{n} X_i^2 - 2 \overline{X}_n\frac{1}{n}\sum_{i=1}^{n} X_i + \overline{X}_n^2  \\
 &= \frac{1}{n}\sum_{i=1}^{n} X_i^2 - 2 \overline{X}_n^2  + \overline{X}_n^2  \\
 &= \frac{1}{n}\sum_{i=1}^{n} X_i^2 - \overline{X}_n^2.
 \end{align*}
$$

The second term in the expression for $S_n^2$ is
determined by

$$
\begin{align*}
    \overline{X}_n&\xrightarrow[n \to \infty]{\mathbb{P}}\mathbb{E}\left[X\right] &\quad&\text{(LLN)} \\
    \mathbb{E}\left[X\right] &= \mu &\quad&\text{(given)}. \\
    g(\overline{X}_n) &\xrightarrow[n \to \infty]{\mathbb{P}}g(\mu) &\quad&\text{(CMT)} \\
    \implies \overline{X}_n^2 &\xrightarrow[n \to \infty]{\mathbb{P}}\mu^2.
\end{align*}
$$

The first term in the expression for $S_n^2$ is then
determined by

$$
\begin{align*}
    \frac{1}{n}\sum_{i=1}^{n} X_i^2 &\xrightarrow[n \to \infty]{\mathbb{P}}\mathbb{E}\left[X^2\right] &\quad&\text{(LLN)} \\
    \operatorname{Var}\left(X\right) &= \mathbb{E}\left[X^2\right] - \mathbb{E}\left[X\right]^2 &\quad&\text{(definition)} \\
    \implies \mathbb{E}\left[X^2\right] &= \operatorname{Var}\left(X\right) + \mathbb{E}\left[X\right]^2 \\
                     &= \sigma^2 + \mu^2. &\quad&\text{(given)} \\
    \therefore S_n^2 &\xrightarrow[n \to \infty]{\mathbb{P}}\sigma^2 + \mu^2 - \mu^2 \\
    \implies S_n^2 &\xrightarrow[n \to \infty]{\mathbb{P}}\sigma^2 \tag*{◻}
  \end{align*}
$$

</div>

Thus,
$V \xrightarrow[n \to \infty]{\mathbb{P}}\frac{n \sigma^2}{\sigma^2} = n$,
a constant.

<div class="theorem" markdown=1>

**Theorem 6** (Slutsky’s Theorem). *If the sequences of random variables
$X_n~\xrightarrow[n \to \infty]{(d)}~X$, and
$Y_n~\xrightarrow[n \to \infty]{(d)}~c$, a constant, then*

$$
\begin{align*}
    X_n + Y_n &\xrightarrow[n \to \infty]{(d)}X + c \text{, and} \\
    X_n Y_n &\xrightarrow[n \to \infty]{(d)}cX.
  \end{align*}
$$

</div>

Since convergence in probability implies convergence in distribution,
and $Z
\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right)$,
Slutsky’s theorem implies that

$$
\begin{align*}
  t_{n-1} = \frac{\displaystyle Z}{\displaystyle \sqrt{\frac{V}{n-1}}} 
  &\xrightarrow[n \to \infty]{(d)}\frac{\displaystyle \mathcal{N}\left( 0, 1 \right)}{\displaystyle \sqrt{\frac{n}{n-1}}} \\ \implies t_{n-1} 
  &\xrightarrow[n \to \infty]{(d)}\mathcal{N}\left( 0, 1 \right).  \tag*{■}
\end{align*}
$$

</div>
</div>

## Conclusions
Since we've failed to reject our null hypothesis that being happy and being in
a relationship are independent, we might be tempted to cite the study as proof
that you do not need to be in a relationship to be happy. Let's not forget
though that by failing to reject the null hypothesis, we haven't actually
*proven* anything, we've just *failed to prove* otherwise. How might we redesign
our study if we wanted to have a stronger result?

The sample size of $N = 1000$ is fairly large, but what might be some other
confounding factors that our study failed to capture? For one, we do not know
anything about the homogeneity of the population from which the sample was
taken. This data could have been taken strictly from a culture that values
individuality over companionship, in which case it is difficult to extrapolate
any conclusions to humanity as a whole. 

Looks like we need to collect more evidence. Until then, I'd say it's safer for
our mental health to assume that happiness and being in a relationship are
indeed independent!


[ocw]: https://ocw.mit.edu/courses/mathematics/18-650-statistics-for-applications-fall-2016/index.htm
