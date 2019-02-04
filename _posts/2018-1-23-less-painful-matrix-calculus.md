---
layout: post
title: Less Painful Matrix Calculus
use_math: true
---

<div style="display none">
$$
  \newcommand{\bm}[1]{\boldsymbol{#1}}
  \newcommand{\mb}[1]{\mathbb{#1}}
  \newcommand{\mc}[1]{\mathcal{#1}}
  \newcommand{\trace}{\text{trace}}
  \newcommand{\real}{\mathbb{R}}
  \renewcommand{\vec}{\text{vec}}
  \newcommand{\norm}[1]{\left\Vert #1 \right\Vert}
$$
<div/>

---
**Quick warning:** I'm writing this post mostly as a reference since I've had a hard time finding a good at-a-glance resource for differentials. 
In other words, no apologies for dryness.

If you already have some familiarity with this stuff and are just trying to refresh, feel free to skip to the bottom for examples.

---

There have been various points throughout my grad school experience where I've found myself struggling to compute derivatives of some fancy-schmancy matrix functions like nuclear norms $\norm{ \bm{A} }\_\*$
or even atomic norms $\norm{\bm{A}}\_{\mc{A}}$ (which is a topic for another blog post).
These kinds of problems used to induce a lot stress googling until about a year ago when *my whole life changed*. 
Oddly reminiscent of the first time I used a bidet.
Well, not really. 
But at the very least, I got to spend less time reading [math overflow posts](https://mathoverflow.net) and more time stuck on the big problems like "what am I making for dinner tonight?"

The method to which I'm referring is differentials. I'm sure that this stuff is fairly well-known, I just have a hard time finding satisfactory references ergo this blog post.
You can find this stuff in [Matrix Calculus](https://www.cambridge.org/core/books/matrix-algebra/BCE8FD2D62006D4061F88E02615B5622) by Karim M. Abadir -- I'll be using some examples found there as well as a set of lecture slides by Gonggou Tang for his convex optimization class (I'm not sure if these slides are publicly available). 
Also, [this paper](https://tminka.github.io/papers/matrix/) by Thomas Minka is great and also covers a lot of this stuff -- check there if you can't find something you want here.

## Definitions

Differentials aren't really anything fancy. 
The concept pretty much comes straight out of the first principles definition of a derivative for a univariate function $f : \real \rightarrow \real$

$$
  f'(x) := \lim\limits_{dx \rightarrow \infty} \frac{f(x + dx) - f(x)}{dx}.
$$

By rearranging, we get

$$
  f(x + dx) = f(x) + f'(x) (dx) + r_x(dx),
$$

where $r_x(dx) - \rightarrow 0$ as $dx \rightarrow 0$. 
Okay, but what does that get us? 
Well, as it turns out, the **differential** of $f$ is defined as

$$
  df(x ; dx) := f'(x) (dx).
$$

This means that if we can easily compute the differential, then we get the derivative for free. 
More generally, for a vector function $f : \real^n \rightarrow \real^m$, we write the differential as 

\begin{equation}\label{eq:1st differential}
  df(x ; dx) := \bm{D}(x) (dx)
\end{equation}

where 

$$
  \bm{D}(x) = \begin{bmatrix}
    \frac{\partial f_1}{\partial x_1} & \cdots & \frac{\partial f_1}{\partial x_n} \\
    \vdots & \ddots & \vdots \\
    \frac{\partial f_m}{\partial x_1} & \cdots & \frac{\partial f_m}{\partial x_n}
  \end{bmatrix} \in \real^{m \times n}
$$

is the derivative matrix. 
The "$;$" notation becomes a bit cumbersome in places, so I'll only use it for definitions and use $df(x)$ otherwise.

As it turns out, differentials work similarly for matrix valued functions $f : \real^{m \times n} \to \real^{l \times k}$, but I'll spare you the details. 
It's worth pointing out that $\bm{D}(x) = \bm{J}(x)^\top$ -- the Jacobian matrix. 
In the scalar-valued case $f : \real^m \to \real$, this means that the derivative vector is simply the gradient transposed $\bm{D}(x) = \nabla f(x)^\top$.

### Second Derivatives

Differentials become particularly helpful once we realize that we can compute second order derivatives just as easily.
Recall the Taylor series expansion

$$
  f(x + dx) = f(x) + \bm{D}(x)(dx) + \frac{1}{2}(dx)^\top \bm{H}(x) (dx) + \mc{O}(\norm{ dx }_2^2),
$$

which gives us a matrix encapsulating second order derivative information $\bm{H}(x)$ called the Hessian matrix.

Given how intimiately connected the Taylor series and differentials are, we should expect there to be a way to get $\bm{H}(x)$ that's just as easy as getting $\bm{D}(x)$.
As it turns out, all we have to do is compute the differential of the differential and match up terms:

\begin{equation}\label{eq:2nd differential}
  d^2 f(x; dx) := (dx)^\top \bm{B}(x) (dx) \quad\Rightarrow\quad \bm{H}(x) = \frac{1}{2}(\bm{B}(x) + \bm{B}(x)^\top).
\end{equation}

Notice that the second differential gives us a matrix $\bm{B}(x)$ -- this matrix doesn't have to be symmetric in general.
This is a problem since the Hessian $\bm{H}(x)$ is by definition symmetric.
This is easily solved by taking the symmetric part of $\bm{B}(x)$, which is precisely what we have in Eq. (\ref{eq:2nd differential}).

## Rules

The above definitions are all well and good, but they don't give us a way to actually find differentials. 
So far, there's no perceivable benefit to using differentials at all! 
Well, just like in kindergarten calculus, we have a standard set of rules for differentials -- the product rule, chain rule, quotient rule, and so forth.
These rules work pretty much exactly the same as in regular calculus, which is great news because it makes the overhead of using them minimal.

### Chain Rule

Let's get the chain rule out of the way first. Remember that the standard vector calculus chain rule is

$$
  \bm{D} f(g(x)) = \bm{D}_f \, f(g(x)) \bm{D}_g \, g(x).
$$

The chain rule for differentials is analogous, it just has a fancy syntax and fancier name -- **Cauchy's Rule of Invariance**:

$$
  d h(x; dx) = df\left( g(x); dg(x; dx) \right), \quad\quad h(x) := f(g(x)).
$$

This looks mysterious, but it's functionally the same thing.
For example, if $f(x) = e^{x^2}$,

$$
  df(x) = d\left( e^{x^2} \right) = e^{x^2} (d(x^2)) = e^{x^2} (2x \cdot dx),
$$

which is exactly what we get when computing derivatives the normal way (differentials aren't too special in 1D).
So that works.
What doesn't work is trying to use Cauchy's invariance for the second differential $d^2$. 
This should be fairly obvious, but heres an example using the same function: if we naively use the chain rule on $d^2$, we get

$$
  d^2\left( e^{x^2} \right) = e^{x^2} (d^2(x^2)) = e^{x^2} (4x^2) (dx)^2,
$$

which is incorrect. In reality,

$$
\begin{aligned}
  d^2 f(x) &= d\left( d\left( e^{x^2} \right) \right) \\
    &= d\left( 2x e^{x^2} \cdot dx \right) \\
    &= 2(dx) e^{x^2} (dx) + 2x e^{x^2} (2x \cdot dx) (dx) \\
    &= \left( 2e^{x^2} + 2x e^{x^2} \right) (dx)^2
\end{aligned}
$$

Note that I've used the product rule here without actually defining it for you.
This is okay since it works exactly the same as in standard calculus.

The important point here is that if you want to compute second order derivatives with differentials, it's easiest to compute the differential of the first differential $d( df(x) )$.
Apparently there are rules set up for second differentials explicitly, but I've never looked them up and I probably never will.

### The Other Rules and Identification Tables

I've mentioned this a couple times, but the standard calculus rules translate pretty much directly over to differentials. 
Here's a summary of the rules (Figures 1, 2) that I'm ripping directly out of Gongguo Tang's lecture slides.

[//]: <> ![table 1]({{ site.url }}/images/1-23-19/rules_table1.PNG)
![](https://github.com/jhelland/jhelland.github.io/blob/master/images/1-23-19/rules_table1.PNG?raw=true)
<center>Figure 1: linearity, product, quotient, and transpose rules</center>

[//]: <> ![table 2]({{ site.url }}/images/1-23-19/matrix_rules.PNG)
![](https://github.com/jhelland/jhelland.github.io/blob/master/images/1-23-19/matrix_rules.PNG?raw=true)
<center>Figure 2: matrix differential rules</center>

In Figure 2, $\otimes$ means the Kronecker product and $\odot$ means the Hadamard product.
Although the chain rule isn't listed in these figures, you can just refer to the previous section.
These rules wouldn't be very useful without identification tables, however, so Figures 3 and 4 have those.
I haven't been able to find a good succinct source of tables like this anywhere, so I thought it would be useful to have those handy.

[//]: <> ![table 3]({{ site.url }}/images/1-23-19/first_identification_table.PNG)
![](https://github.com/jhelland/jhelland.github.io/blob/master/images/1-23-19/first_identification_table.PNG?raw=true)
<center>Figure 3: the first identification table</center>

[//]: <> ![table 4]({{ site.url }}/images/1-23-19/second_identification_table.PNG)
![](https://github.com/jhelland/jhelland.github.io/blob/master/images/1-23-19/second_identification_table.PNG?raw=true)
<center>Figure 4: the second identification table</center>

In Figure 4, we've got this weird subscript $v$ notation. That just means a sort of blockwise transpose:

$$
  \bm{B} = \begin{bmatrix}
    \bm{B}_1 \\ \vdots \\ \bm{B}_m
  \end{bmatrix}, \quad\quad (\bm{B}^\top)_v = \begin{bmatrix}
    \bm{B}_1^\top \\ \vdots  \\ \bm{B}_m^\top
  \end{bmatrix}.
$$

It's usually easy enough to remember the identifications for the first differentials, but the second differentials can be trickier to remember -- especially for matrix-valued functions.

For matrix-valued functions, you probably notice some heavy usage of the $\vec(\cdot)$ and $\otimes$ operators. These have some properties that are useful to know when working in this context (a by no means complete list):

$$
  \trace\left(\bm{A}^\top \bm{B}\right) = \vec(\bm{A})^\top \vec(\bm{B})
$$

for any compatibly sized matrices $\bm{A}$ and $\bm{B}$. 
In addition, let $\bm{A}, \bm{B}, \bm{C}$ be $k \times l$, $l \times m$, and $m \times n$ respectively. 
Then

$$
\begin{aligned}
  \vec(\bm{ABC}) &= (\bm{C}^\top \otimes \bm{A}) \vec(\bm{B}) \\[1em]
  \vec(\bm{ABC}) &= (\bm{I}_n \otimes \bm{AB}) \vec(\bm{C}) = (\bm{C}^\top \bm{B}^\top \otimes \bm{I}_k) \vec(\bm{A}) \\[1em]
  \vec(\bm{AB}) &= (\bm{I}_m \otimes \bm{A}) \vec(\bm{B}) = (\bm{B}^\top \otimes \bm{I}_k)\vec(\bm{A}).
\end{aligned}
$$

Finally, one of the most useful tricks to know:

$$
  \trace( \bm{A} \bm{B} ) = \trace(\bm{B} \bm{A})
$$

assuming that the dimensions line up appropriately. 
There's other properties of the trace like invariance to transposes that are nice to have on hand, but these properties will get you pretty far with differentials.

## Examples

What's a good math reference without some examples?
When I inevitably forget this stuff again, examples will be the first thing I look for.
Granted, there's already some brief examples above, but those are in 1D and don't do much to show off the power of differentials.

### Problem 1

Consider 

$$
  f(x) = \frac{1}{2} x^\top \bm{A} x, \quad\quad \bm{A} \in \real^{n \times n}, \ x \in \real^n
$$

where $\bm{A}$ is symmetric and suppose that we want to compute the gradient and Hessian. Then

$$
  df(x) = \frac{1}{2}( (dx)^\top \bm{A} x + x^\top \bm{A} \cdot dx ) = x^\top \bm{A} \cdot dx.
$$

Using the first identification table, we get $\nabla f(x) = \bm{A}^\top x = \bm{A} x$. 

To get the Hessian, simply apply another differential:

$$
  d^2 f(x) = (dx)^\top \bm{A} \cdot dx,
$$

which gives the Hessian $\bm{H} f(x) = \nabla^2 f(x) = \frac{1}{2}(\bm{A} + \bm{A}^\top) = \bm{A}$. That was easy enough.

### Problem 2

Now suppose we want the gradient and Hessian of

$$
  f(x) = \frac{1}{2} \norm{y - \bm{A}x}_2^2, \quad\quad y \in \real^m, \ \bm{A} \in \real^{m \times n}, \ x \in \real^n.
$$

We could notice that $\norm{ \cdot }_2^2$ is simply an inner product and expand accordingly, but that seems unnecessary. 
Instead, let's use the chain rule.

$$
\begin{aligned}
  df(x) &= \frac{1}{2} d\left( \norm{y - \bm{A}x}_2^2 \right) (dy - d(\bm{A}x)) \\[.5em]
    &= \frac{1}{2} \cdot 2(y - \bm{A}x)^\top(-\bm{A} \cdot dx) \\[.5em]
    &= (\bm{A}x - y)^\top \bm{A} \cdot dx
\end{aligned}
$$

since $d(x^\top x) = 2x^\top \cdot dx$. The gradient is then $\nabla f(x) = \bm{A}^\top (\bm{A}x - y)$.
The second differential is then

$$
  d^2 f(x) = d\left( x^\top \bm{A}^\top \bm{A} - y^\top \bm{A} \right) = (dx)^\top \bm{A}^\top \bm{A} \cdot dx,
$$

which means that $\nabla^2 f(x) = \bm{A}^\top \bm{A}$.

---

More examples coming soon.

<!--
  ### Problem 3

  Alright, let's get funkier. Take

  $$
    f(\bm{U}, \bm{V}) = \frac{1}{2} \norm{ \bm{Y} - \bm{UV}^\top }_F^2, \quad\quad \bm{Y} \in \real^{m \times n}, \ \bm{U} \in \real^{m \times r}, \ \bm{V} \in \real^{n \times r}.
  $$

  where $r > 1$.
  In the previous problems, it would have been simple to compute entrywise derivatives and do some bookkeeping. Here, the bookkeeping starts to become more of a headache. 
  Instead of that nonsense, we can do differentials as follows:

  $$
  \begin{aligned}
    d_{\bm{U}} f(\bm{U}, \bm{V}) &= \frac{1}{2} \trace\left( (-d\bm{U} \cdot \bm{V}^\top)(\bm{UV}^\top - \bm{Y})^\top + (\bm{Y} - \bm{UV}^\top)(-d\bm{U} \cdot \bm{V}^\top)^\top \right) \\[.5em]
      &= \trace\left( d\bm{U} \cdot \bm{V}^\top(\bm{UV}^\top - \bm{Y})^\top \right) \\[.5em]
      &= \vec((\bm{UV}^\top - \bm{Y})\bm{V})^\top \vec(d\bm{U}),
  \end{aligned}
  $$

  which by identification using Figure 3 gives $\nabla_{\bm{U}} f(\bm{U}, \bm{V}) = \vec\left( (\bm{UV}^\top - \bm{Y}) \bm{V} \right)$. Note that by $\nabla_{\bm{U}}$, I mean the gradient with respect to $\bm{U}$ only, in case that wasn't clear. 

  You can do some similar algebra and find that $\nabla_{\bm{V}} f(\bm{U}, \bm{V}) = \vec\left( (\bm{UV}^\top - \bm{Y})^\top \bm{U} \right)$. This means that the overall gradient is

  $$
    \nabla f(\bm{U}, \bm{V}) = \begin{bmatrix}
      \vec\left( (\bm{UV}^\top - \bm{Y}) \bm{V} \right) \\[.5em]
      \vec\left( (\bm{UV}^\top - \bm{Y})^\top \bm{U} \right)
    \end{bmatrix} \in \real^{(m + n)r}.
  $$

  Having not done this problem the hard way (and having no immediate plans to), I'm going to go out on a limb and say that this approach is much quicker than trying to compute derivatives of individual components.

  "But wait!" cries the interested reader. We haven't computed the Hessian yet! Fear not, dear reader, for that is precisely what we'll do next.
  In the same way that we computed blocks of $\nabla f$ and then combined, we'll need to compute four blocks of the Hessian $\nabla^2 f$ before we have the final answer. 
  We'll start as before: with $d_{\bm{U}}$.

  $$
  \begin{aligned}
    d_{\bm{U}}^2 f(\bm{U}, \bm{V}) &= \vec\left( d\bm{U} \cdot \bm{V}^\top \bm{V} \right)^\top \vec( d\bm{U} ) \\
      &= \vec(d\bm{U})^\top \left( \bm{V}^\top \bm{V} \otimes \bm{I}_m \right) \vec(d \bm{U}),
  \end{aligned}
  $$

  so we get $\nabla_{\bm{U}}^2 f(\bm{U}, \bm{V}) = \bm{H}_{\bm{UU}} = \left(\bm{V}^\top \bm{V} \otimes \bm{I}_m \right) \in \real^{mr \times mr}$.
  I'm tracking dimensions a little more carefully now because once Kronecker products get involved, it's easy to lose track of what's happening.
  Additionally, simply matching up dimensions is my favorite way of double checking my linear algebra work in general.

  Okay, so we can do similar algebra and get the other main diagonal block as $\bm{H}_{\bm{VV}} = \left( \bm{U}^\top \bm{U} \otimes \bm{I}_n \right) \in \real^{nr \times nr}$. 
  The only thing left is to compute the off-diagonal blocks. Fortunately, we can recall that the Hessian kkk
-->
