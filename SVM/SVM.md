## 1.拉格朗日对偶性

&emsp;&emsp;在约束最优化问题中，常常利用拉格朗日对偶性（Lagrange duality）将原始问题转换为对偶问题，通过求解对偶问题而得到原始问题的解。

### 1.1 原始问题

&emsp;&emsp;假设 $f(x),c_i(x),h_j(x)$ 是定义在 $\textbf{R}^n$ 上的连续可微函数。考虑约束最优化问题：
$$
min_{x\in \textbf{R}^n}f(x)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.1  \\ s.t.\ \ c_i(x)\leq0,\ \ \ i=1,2,···,k \ \ \ \ \ \ \ \ \ \ 1.2\\ \ \ \ \ \ \ \ \ h_j(x)=0,\ \ \ j=1,2,···,l\ \ \ \ \ \ \ \ \ \ \ 1.3
$$
将这个约束最优化问题称为原始最优化问题或是原始问题。

&emsp;&emsp;首先，引进**广义拉格朗日函数**（generalized Langrange function）:
$$
L(x,\alpha,\beta)=f(x)+\sum_{i=1}^{k}\alpha_ic_i(x)+\sum_{j=1}^{l}\beta_jh_j(x)\ \ \ \ \ \ \ \ \ \ 1.4
$$
其中，$x=(x^{(1)},x^{(2)},···,`x^{(n)})^T\in\textbf{R}^n$，$\alpha_i,\beta_j$ 是拉格朗日乘子，其中，$\alpha_i\geq0$。这里我们首先考虑 $x$ 的函数：
$$
\theta_p(x)=max_{\alpha,\beta:\alpha_i\geq 0}L(x,\alpha,\beta)\ \ \ \ \ \ \ \ \ \ 1.5
$$
上面公式中的下标 $p$ 指的是原始问题。

&emsp;&emsp;假设给定某个 $x$。如果 $x$ 违反原始问题的约束条件，也就是存在某个 $i$ 使得 $c_i(w)>0$ 或者存在某个 $j$ 使得 $h_j(x)\neq0$，那么则有：
$$
\theta_p(x)=max_{\alpha,\beta:\alpha_i\geq 0}\left[f(x)+\sum_{i=1}^{k}\alpha_ic_i(x)+\sum_{j=1}^{l}\beta_jh_j(x)\right]=+\infty\ \ \ \ \ \ \ \ \ \ 1.6
$$
因为如果某个 $i$ 使得约束 $c_i(x)>0$，则可令 $\alpha_i\rightarrow+\infty$，如果某个 $j$ 使得 $h_j(x)\neq0$，则可令 $\beta_j$ 使 $\beta_jh_j(x)\rightarrow+\infty$，而将其余各 $\alpha_i,\beta_j$ 均取值为 0。

&emsp;&emsp;相反的，如果 $x$ 满足约束条件（1.2）和（1.3），则由（1.5）和（1.4）可知：$\theta_p(x)=f(x)$，因此可以得到：
$$
\theta_p(x)=\begin{equation}  
\left\{  
             \begin{array}{lr}  
             f(x) ,\ \ \ \ \ \ \ \ x满足原始问题约束\\  
            +\infty\ \ \ \ \ \ \ \ \ \   其他\\  
             \end{array}  
\right.  
\end{equation} 
$$
所以如果考虑极小化问题：
$$
min_x\theta_p(x)=min_x\ \ max_{\alpha,\beta:\alpha_i\geq 0}L(x,\alpha,\beta)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.8
$$
上面的等式表明它与原始最优化问题（1.1）~（1.3）是等价的，也就是它们由相同的解。问题$min_x\ \ max_{\alpha,\beta:\alpha_i\geq 0}L(x,\alpha,\beta)$ 称为 **广义拉格朗日的极小极大** 问题。这样一来，就把原始最优化问题表示为广义拉格朗日函数的极小极大问题。为了接下来的推导方便，这里定义原始问题的最优值：
$$
p^*=min_x\theta_p(x)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.9
$$
称为原问题的值。



### 1.3 对偶问题

&emsp;&emsp;定义：
$$
\theta_D(\alpha,\beta)=min_xL(x,\alpha,\beta)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.10
$$
再考虑极大化 $\theta_D(\alpha,\beta)=min_xL(x,\alpha,\beta)$，也就是：
$$
max_{\alpha,\beta:\alpha_i\geq 0}\theta_D(\alpha,\beta)= max_{\alpha,\beta:\alpha_i\geq 0}\ \ min_xL(x,\alpha,\beta)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.11
$$
上面的问题称为广义拉格朗日函数的极大极小问题。

&emsp;&emsp;可以将广义拉格朗日函数的极大极小问题表示为约束最优化问题：
$$
max_{\alpha,\beta:\alpha_i\geq 0}\theta_D(\alpha,\beta)=max_{\alpha,\beta}min_xL(x,\alpha,\beta)\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ 1.12\\s.t.\ \ \ \ \ \ \ \alpha_i\geq0,\ \ \ \ i=1,2,···,k\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ 1.13
$$
称为原始问题的对偶问题，定义对偶问题的最优值：
$$
d^*=max_{\alpha,\beta:\alpha_i\geq0}\theta_D(\alpha,\beta)\ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ 1.14
$$
称为对偶问题的值。



### 1.4 原始问题和对偶问题的关系

&emsp;&emsp;这里主要介绍几条重要的定理：

**定理1**

> 若原始问题与对偶问题都有最优值，则有：
> $$
> d^*=max_{\alpha,\beta:\alpha_i\geq0}min_xL(x,\alpha,\beta)\leq min_x\ \ max_{\alpha,\beta:\alpha_i\geq 0}L(x,\alpha,\beta)=p^*\ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ \ \ \ 1.15
> $$
>



**推论1**

> 设 $x^*$ 和 $\alpha^*,\beta^*$ 分别是原始问题（1.1）~（1.3）和对偶问题 （1.12）~（1.13）的可行解，且有**$d^*=p^*$**，则 $x^*$ 和 $\alpha^*,\beta^*$ 分别是原问题和对偶问题的**最优解**。



**定理2**

> 考虑原始问题（1.1）~（1.3）和对偶问题 （1.12）~（1.13），假设函数 $f(x),c_i(x)$ 都是凸函数，$h_j(x)$ 是仿射函数；并且假设不等式约束 $c_i(x)$ 是严格可行的，也就是不存在等号的情况————
>
> - 则存在 $x^*$ 和 $\alpha^*,\beta^*$ ，使得 $x^*$ 是原问题的解，$\alpha^*,\beta^*$ 是对偶问题的解，并且有：
>
> $$
> p^*=d^*=L(x^*,\alpha^*,\beta^*)
> $$
>



**定理3**

