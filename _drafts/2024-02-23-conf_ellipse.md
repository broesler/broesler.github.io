---
layout: post
title:  "kD Gaussian Confidence Ellipses"
date:   "2024-02-23 21:10 -0500"
categories: statistics
tags: statistics python
reading_time: 15
excerpt: What is the confidence region for a kD Gaussian distribution? And how
    do we plot it? 
---

<!-- MARGIN NOTES -->
<!-- <span --> 
<!-- class="marginnote-wrapper"><label for="sn-0" class="margin-toggle">&#8862;</label> -->
<!-- <input type="checkbox" id="sn-0" class="margin-toggle"/> -->
<!-- <span class="marginnote"> -->
<!-- Note that the two samples may have different sizes (if -->
<!-- $n \ne m$). -->
<!-- </span></span> -->

I recently went to make a simple plot of the confidence ellipses for a 2D
Gaussian distribution fit to some data. I had computed contours of the pdf
numerically, but then found this [matplotlib example][matplotlib-example] which
linked to a nice [blog post][carsten-blog] describing the basic parameters of
a confidence ellipse. The problem was, I couldn't get the two solutions to
agree. What follows is my investigation into the "bug", and a richer depth of
material on multivariate normal distributions than I could have imagined.

## Trouble Brewing

I discovered this problem when I plotted a small sample from simulated
data, as shown in
<a href="#fig:initial" data-reference-type="ref" data-reference="fig:initial">Figure 1</a>.

<figure>
<img src="{{ '/assets/images/conf_ellipse/initial_problem_notes.pdf' | absolute_url }}"
id="fig:initial"  /><figcaption><span class="fig_number">Figure 1</span>. 
The initial small dataset with a "1$\sigma$" confidence ellipse plotted
numerically, and the supposed "1$\sigma$" confidence ellipse using the
``confidence_ellipse`` function described in Carsten's blog and the
corresponding matplotlib example. Points outside of the numerical ellipse are
colored red.
</figcaption>
</figure>

The trouble was, who was wrong? One the one hand, Carsten's solution had made it
into the matplotlib examples. Surely there are some maintainers with a wealth of
statistical knowledge who would have vetted these examples, or at least some
dedicated citizens of the internet. On the other hand, my numerical method only
occupied four lines of code (five if you count the contour label!). It does
rely upon the matplotlib `contour` function, but short of delving into that
machinery, I couldn't see where I had gone astray.

## Covariance Confusion

Before getting into the math notation, let's clear up some terminology.

Carsten's blog uses the terms "covariance" ellipse and "confidence" ellipse
somewhat interchangeably, other [blog posts][error-ellipse-blog] use terms like
"error ellipse, a.k.a. confidence ellipse", and other terms such as "standard
deviational" ellipse or "tolerance region" appear in academic
[literature][wang-SDHE].

These terms all refer to contours of constant density of the multivariate normal
distribution, which turn out to be ellipses that are centered on the mean. The
``confidence_ellipse`` function on the matplotlib example (as described in
Carsten's blog) allows you to set the scale of the ellipse by the argument
``n_std`` for "number of standard deviations". We would alternatively like to be
able to specify the scale of the ellipse by its "level" or the total probability
mass that it encloses.

You are probably familiar with the 1-D version of this idea, known as the
["3$\sigma$" rule][3-sigma], which states that the probability that a random
value drawn from a normal distribution lies within 1, 2, or 3 standard
deviations of the mean is 68.3%, 95.4%, and 99.7%, respectively. Given a mean
value, you can specify an interval about the mean by either the number of
standard deviations, or by the percentage level it encloses.

Plenty has been written about the many
[misinterpretations of confidence intervals][misinterp-journal], so I will offer
the following disclaimer: Confidence intervals (or *regions* in multiple
dimensions) are *random* intervals constructed by a statistical procedure with
respect to a particular sample taken in an experiment for purposes of
statistical inference. The ellipses described in
the remainder of this blog post should be
only be interpreted as "confidence regions" at your own discretion in the
context of your specific data and experimental design. 

## Problem Statement

We can now more precisely define the problem we are trying to solve.

Let $\mathbf{X} = \[X_1, \dots, X_k\]^\T$ be a random vector drawn i.i.d. from
a multivariate normal distribution 

$$
\mathbf{X} \sim \mathcal{N}(\boldsymbol{\mu}, \Sigma)
$$

where $\boldsymbol{\mu}$ is the $k$-dimensional mean vector 

$$
\boldsymbol{\mu} = \E{\mathbf{X}},
$$

and $\Sigma$ is the $k \times k$ covariance matrix

$$
\Sigma_{i,j} = \Cov{X_i, X_j} = \E{(X_i - \mu_i)(X_j - \mu_j)}
$$

for $i, j \in \\{1, 2, \dots, k\\}$.

We want to define a region $\mathcal{R}$ such that 
$\Prob{\mathbf{X} \in \mathcal{R}} \ge \alpha$ for some $0 \le \alpha \le 1$.
Define this region as

$$
\mathcal{R} \equiv \left\{ \mathbf{X} : R(\mathbf{X}; \boldsymbol{\mu}, \Sigma)^2 \le C^2 \right\}
$$

for some function $R$ that represents a distance from the mean, and constant
$C$. Squaring these values makes the inequality easier to work with since they
are always positive. We will now set out to define these values and then we will
have our desired region.

## The Distance Function

In 1-D, the concept of distance from the mean in terms of "number of standard
deviations" is embodied in the $z$-value. Given a normal random value 

$$
X \sim \mathcal{N}\left(\mu, \sigma^2\right),
$$

we center and normalize it to get a standard normal variable

$$
Z = \frac{X - \mu}{\sigma} \sim \mathcal{N}(0, 1).
$$

How should we extend this idea into multiple dimensions? Should we standardize
each dimension individually? It turns out the $k$-D analogy of this
normalization is

$$
R^2 = (\mathbf{X} - \boldsymbol{\mu})^\T \Sigma^{-1} (\mathbf{X} - \boldsymbol{\mu}).
$$

The value $R$ is known as the *Mahalanobis distance*.<span
class="sidenote-wrapper"><label for="sn-1" class="margin-toggle sidenote-number"></label><input type="checkbox" id="sn-1" class="margin-toggle"/><span
class="sidenote">
If our data is real-valued and has non-zero variance
in all dimensions, the covariance matrix will be positive definite, so its
inverse exists and is also positive definite.
</span></span>
This standardization using the entire covariance matrix retains all of the
necessary information about the correlation between each dimension. If we were
to standardize each dimension individually, we would only be using the diagonal
of $\Sigma$, and throwing away all of that additional information about the
relationship between each dimension.

We will use this distance function to define our region:

$$
\mathcal{R} \equiv \left\{ \mathbf{X} :
(\mathbf{X} - \boldsymbol{\mu})^\T \Sigma^{-1} (\mathbf{X} - \boldsymbol{\mu})
\le C^2 \right\}.
$$

## Computing Probabilities

Now that we have a distance function, we need to determine the maximum squared
distance $C^2$ so that the region encloses the desired amount of probability.
Notice that this distance $R^2$ we have defined is
itself a random variable, so we can ask the question: how is it distributed? 

Combining the 1-D analogy with some vector multiplication rules, we can see that 
$R^2$ is actually a *sum* of $k$ squared, independent, standard normal
random variables, which we know to be distributed according to the chi-squared
distribution with $k$ degrees of freedom. Thus, the probability that a random
vector $\mathbf{X}$ is within a distance $r$ of the mean is

$$
\Prob{R(\mathbf{X}) \le r} = \Prob{R(\mathbf{X})^2 \le r^2} = F_{\chi^2_k}(r^2)
$$

where $F_{\chi^2_k}(r^2)$ is the cumulative distribution function (CDF) of the
chi-squared distribution. Since we want to specify this probability, we can
invert the CDF

$$
\begin{align*}
\Prob{R(\mathbf{X})^2 \le r^2} = 1 - \alpha &= F_{\chi^2_k}(r^2) \\
\chi^2_k(1 - \alpha) &= r^2,
\end{align*}
$$

where $\chi^2_k(q)$ is the inverse of the CDF, the *quantile function*, for
probability $q$ of the chi-squared distribution.
We can thus let our constant $C^2$ be this value, giving our definition
of the region

$$
\mathcal{R} \equiv \left\{ \mathbf{X} :
(\mathbf{X} - \boldsymbol{\mu})^\T \Sigma^{-1} (\mathbf{X} - \boldsymbol{\mu})
\le \chi^2_k(q) \right\}
$$

for level $q = 1 - \alpha$.

Note that we have defined this region assuming we know the population mean and
covariance matrix. If either the mean or covariance are unknown, we can use the
typical unbiased estimates:

$$
\begin{align*}
\boldsymbol{\mu} &\approx \mathbf{\bar{X}} = \frac{1}{N} \sum_{i=1}^N \mathbf{X}_i \\
\Sigma &\approx \hat{\Sigma} = \frac{1}{N-1} \sum_{i=1}^N \left( \mathbf{X} - \mathbf{\bar{X}} \right)
    \left( \mathbf{X} - \mathbf{\bar{X}} \right)^\T
\end{align*}
$$

A more complete proof that this region is a true confidence region (indeed it is
the smallest possible region that encloses probability $1 - \alpha$), as well as
proof that the region defined by using these estimates converges in probability
to the true region is given in [Siotani (1964)][siotani-1964].

## Converting Between Levels and Standard Deviations
$q = \chi^2_k(z^2)$
<a href="#tab:tab1" data-reference-type="ref" data-reference="tab:tab1">Table 1</a>.

<div id="tab:tab1" style="width: auto; margin: auto;" markdown=1>
<span class="marginnote-wrapper"><label for="sn-0" class="margin-toggle">&#8862;</label>
<input type="checkbox" id="sn-0" class="margin-toggle"/>
<span class="marginnote">note here.</span></span>

|   dimensions |$z =$ 1 |      2 |      3 |      4 |      5 |      6 |
|-------------:|-------:|-------:|-------:|-------:|-------:|-------:|
|            1 | 0.6827 | 0.9545 | 0.9973 | 0.9999 | 1.0000 | 1.0000 |
|            2 | 0.3935 | 0.8647 | 0.9889 | 0.9997 | 1.0000 | 1.0000 |
|            3 | 0.1987 | 0.7385 | 0.9707 | 0.9989 | 1.0000 | 1.0000 |
|            4 | 0.0902 | 0.5940 | 0.9389 | 0.9970 | 0.9999 | 1.0000 |
|            5 | 0.0374 | 0.4506 | 0.8909 | 0.9932 | 0.9999 | 1.0000 |
|            6 | 0.0144 | 0.3233 | 0.8264 | 0.9862 | 0.9997 | 1.0000 |
|            7 | 0.0052 | 0.2202 | 0.7473 | 0.9749 | 0.9992 | 1.0000 |
|            8 | 0.0018 | 0.1429 | 0.6577 | 0.9576 | 0.9984 | 1.0000 |
|            9 | 0.0006 | 0.0886 | 0.5627 | 0.9331 | 0.9970 | 1.0000 |
|           10 | 0.0002 | 0.0527 | 0.4679 | 0.9004 | 0.9947 | 0.9999 |

</div>

## The Ellipsoid

What does this region $\mathcal{R}$ actually *look* like? Well, it should be no
secret at this point that it is an ellipsoid (the $k$-D analogue of an ellipse).

We can get some intuition if we look at the probability density function of the
multivariate normal:

$$
f_\mathbf{X}(\mathbf{x}) = \frac{1}{\sqrt{(2\pi)^k |\Sigma|}}
\exp \left(  -\frac{1}{2} (\mathbf{x} - \boldsymbol{\mu})^\T
              \Sigma^{-1} (\mathbf{x} - \boldsymbol{\mu}) \right)
$$

and notice that the argument of the exponential is the squared Mahalanobis
distance (scaled by $-\frac{1}{2}$). Since the outer surface of our region
$\mathcal{R}$ is the set of all vectors where $R^2 = C^2$, a constant, this
means that our region is bounded by an equiprobability surface of the
distribution.

There is a more in-depth explanation of the [geometry of the covariance
matrix][geom-interp-blog].

## Reduction to 2D


## Carsten's Mistake

Carsten's blog post shows an excellent exploration of the relative proportions
and rotation of the 2D ellipse, and he ultimately arrives at the derivation of
the eigendecomposition of the covariance matrix. The mistake he makes is exactly
the standardization issue we discussed above, where he implicitly assumes that
the bounding box of the ellipse should be scaled by the number of standard
deviations of $x$ and $y$ in each dimension separately. As we have discussed, if
you scale the ellipse in this way, it will not enclose the correct amount of the
probability distribution.

I can offer two possibilities as to why Carsten did not catch this issue:

1. The idea of standardizing each dimension individually *does* produce
    a valid cumulative distribution function, in the sense that

    $$
    F(\mathbf{x}) = \Prob{\mathbf{X} \le \mathbf{x}} 
    \text{, where } \mathbf{X} \sim \mathcal{N}\left(\boldsymbol{\mu}, \Sigma\right)
    $$

    describes the probability that all components of $\mathbf{X}$ are less than
    or equal to their corresponding values in $\mathbf{x}$. This definition is
    how Matlab's [``mvncdf``][matlab-mvncdf] function works, and is also
    [reviewed][balakrishnan-2014] in other [literature][genz-2009].
    
    In contrast, we have defined the CDF in terms of the Mahalanobis distance $R$:

    $$
    F(r) = \Prob{R \le r} \text{, where } R^2 \sim \chi^2_k
    $$

    which is also a valid CDF, but correctly scales the ellipse to enclose the
    desired amount of probability.

2. The relative error between the scaling factors computed by these two methods
    is not all that large, especially at greater distances from the mean for
    ellipses that encompass more probability. <a href="#fig:rel_error"
    data-reference-type="ref" data-reference="fig:rel_error">Figure 2</a> shows
    a plot of the relative error between $\sqrt{\chi^2_2(q)}$ and the number of
    standard deviations from the mean, $z$, where $q = 1 - 2\Phi(-z)$, and
    $\Phi(z)$ is the standard normal CDF. Carsten's
    default value is ``n_std`` $= z = 3$. For ellipses beyond about 3$\sigma$ (or "> 95%
    confidence ellipses"), the relative error is below 15%, which wouldn't be
    particularly obvious to the naked eye.
    <a href="#fig:3s" data-reference-type="ref" data-reference="fig:3s">Figure 3</a>
    shows both "3$\sigma$" ellipses, plotted on a sample of $N = 1000$ points
    drawn from the same distribution as in 
    <a href="#fig:initial" data-reference-type="ref" data-reference="fig:initial">Figure 1</a>.
    There are only a handful of points outside of either ellipse at the 
    "99.7% level". The difference would be even more obscured with a sparser
    dataset.

<figure>
<img src="{{ '/assets/images/conf_ellipse/rel_error.pdf' | absolute_url }}"
id="fig:rel_error"  /><figcaption><span class="fig_number">Figure 2</span>. 
Plot of the relative error of two methods to compute the multivariate normal CDF.
</figcaption>
</figure>

<figure>
<img src="{{ '/assets/images/conf_ellipse/three_sigma_wrong_notes.pdf' | absolute_url }}"
id="fig:3s"  /><figcaption><span class="fig_number">Figure 3</span>. 

A larger sample of $N = 1000$ from the same distribution as 
<a href="#fig:initial" data-reference-type="ref" data-reference="fig:initial">Figure 1</a>
with a "3$\sigma$" confidence ellipse plotted
numerically, and the supposed "3$\sigma$" confidence ellipse using the
``confidence_ellipse`` function described in Carsten's blog and the
corresponding matplotlib example. Points outside of the numerical ellipse are
colored red.
</figcaption>
</figure>


To his credit, Carsten does acknowledge this issue in an edit to the post, in
that the "3$\sigma$" rule does not apply exactly in multiple dimensions, but he
falls short in the explanation of why, and in the correction to the problem.
I hope my post has added some clarity to this subtle issue. 


## Conclusions


<p class="message" markdown=1>
The entire source code for the figures and algorithms in this post is available
[here][gaussian_ellipses.py].
</p>


<!-- links -->
[carsten-blog]: https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
[matplotlib-example]: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html
[gaussian_ellipses.py]: https://github.com/broesler/stats_rethinking/blob/master/ch14/gaussian_ellipses.py
[error-ellipse-blog]: https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix
[misinterp-journal]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4877414/
[wang-SDHE]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4358977/
[3-sigma]: https://en.wikipedia.org/wiki/Normal_distribution#Standard_deviation_and_coverage
[siotani-1964]: https://rdcu.be/dzxmD
[geom-interp-blog]: https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
[matlab-mvncdf]: https://www.mathworks.com/help/stats/multivariate-normal-distribution.html
[balakrishnan-2014]: https://doi.org/10.1002/9781118445112.stat01249
[genz-2009]: https://link.springer.com/chapter/10.1007/978-3-642-01689-9_1
