# Early-Time Edge Identifiability and a Later Observation Collision in Hub-Star Kuramoto Networks

**Research No. 83 — public preprint v1.0 (24 November 2026)**  
**Author:** Ryutaro Yonezu  
**Status:** Independent AI mathematical audit PASS (24 November 2026); Windows execution reported and transcript preserved, original native log not supplied; **not human peer reviewed**; novelty **INCONCLUSIVE**; references **PARTIAL**. **All rights reserved.**

## Abstract

We study identification of a designated leaf-to-leaf edge in a hub-star Kuramoto network from one hub-phase observation after a prescribed reset. We derive, directly from the differential equation, a graph-uniform early-time interval on which the hub phases of any two networks with opposite designated-edge states are strictly ordered. We also construct two five-vertex networks with opposite designated-edge states whose hub phases coincide at a later positive time. The collision is established by a directed-rounding MPFR enclosure of explicit-midpoint iterates combined with an analytic global discretization error bound, and was independently checked by an AI auditor using a separate interval implementation and high-accuracy integration. The collision concerns one particular observation time, not all times or multiple-time measurements.

## 1. Model and observables

Fix $N\ge 3$. The vertices are $\{0,1,\ldots,N-1\}$; vertex $0$ is the hub, adjacent to all leaves. Other edges are arbitrary, undirected leaf-to-leaf edges. Let $a_{jk}\in\{0,1\}$ denote leaf adjacency and take $a_{12}$ as the designated unknown edge. We consider the zero-frequency Kuramoto equations

$$
\dot\theta_i(t)=\frac1N\sum_{j\sim i}\sin\bigl(\theta_j(t)-\theta_i(t)\bigr).
$$

We observe the lifted hub phase $h_G(t)=\theta_0(t)$ after the common reset

$$
\theta(0)=(0,\alpha,\pi-\alpha,-\alpha,\ldots,-\alpha),\qquad
\alpha=\arcsin(1/\sqrt3).
$$

The observations are real-valued lifted phases, not wrapped measurements.

## 2. Uniform early-time separation

**Lemma 1 (direct reset identity).** There exists a number $H_N$, depending only on $N$ and the displayed reset, such that for every admissible graph $G$,

$$
N^2h_G''(0)=H_N+\frac{8}{3\sqrt3}\,a_{12}(G).
$$

**Proof.** Differentiating the hub equation once gives

$$
h_G''(0)=\frac1N\sum_{j=1}^{N-1}\cos(p_j-p_0)
\bigl(\dot\theta_j(0)-\dot\theta_0(0)\bigr),
$$

where $p=\theta(0)$ and $p_0=0$. The hub-edge contributions, and the term involving $\dot\theta_0(0)$, depend only on the fixed reset and $N$. For a leaf-to-leaf edge $\{j,k\}$, its contribution to $N^2h_G''(0)$ is

$$
a_{jk}\sin(p_k-p_j)\bigl(\cos p_j-\cos p_k\bigr).
$$

Every pair except $\{1,2\}$ makes this factor vanish: for leaves $j,k\ge3$, both phases are $-\alpha$; for $\{1,k\}$ with $k\ge3$, the cosines agree; and for $\{2,k\}$ with $k\ge3$, the phase difference is $\pi$. For $\{1,2\}$ the factor is

$$
\sin(\pi-2\alpha)\bigl(\cos\alpha-\cos(\pi-\alpha)\bigr)
=2\sin(2\alpha)\cos\alpha=\frac8{3\sqrt3}.
$$

All edge-independent terms can therefore be absorbed into $H_N$. $\square$

**Theorem 1 (unconditional early-time separation).** Let $G_1$ and $G_0$ be any two admissible networks with the same $N$ and $a_{12}(G_1)=1$, $a_{12}(G_0)=0$. Then

$$
h_{G_1}(t)>h_{G_0}(t)\quad\text{whenever}\quad 0<t<\frac1{2\sqrt3\,N^2}.
$$

**Proof.** In the infinity norm, for the vector field $F_G$ and any $\theta$,

$$
\|F_G\|_\infty<1,\qquad \|DF_G\|_\infty<2,\qquad
\|D^2F_G\|_\infty<4.
$$

The Hessian bound follows by estimating each edge contribution using
$|(v_j-v_i)(w_j-w_i)|\le4\|v\|_\infty\|w\|_\infty$.
As $\theta'''=D^2F_G[F_G,F_G]+DF_G\,DF_G\,F_G$, it follows that $|h_G'''(t)|<8$. Write $D=h_{G_1}-h_{G_0}$. The common reset gives $D(0)=D'(0)=0$; Lemma 1 gives $D''(0)=8/(3\sqrt3 N^2)$; and $|D'''(t)|<16$. Taylor's theorem yields

$$
D(t)\ge \frac4{3\sqrt3\,N^2}t^2-\frac83t^3>0
\qquad\left(0<t<\frac1{2\sqrt3\,N^2}\right).
$$

$\square$

This is a separation statement for noise-free means, not a sample-efficiency result. When $t$ is of order $N^{-2}$, the displayed lower bound is of order $N^{-6}$.

## 3. A five-vertex later-time collision

Set $N=5$. Both networks include all hub edges $\{0,i\}$ for $i=1,2,3,4$. The leaf edges of $G_0$ are $\{1,3\},\{1,4\},\{2,3\},\{2,4\}$; the only leaf edge of $G_1$ is $\{1,2\}$. With the common reset above, let $D(t)=h_{G_1}(t)-h_{G_0}(t)$.

**Theorem 2 (computer-assisted collision).** For this graph pair,

$$
D(3.10)>0>D(3.12).
$$

Consequently, there exists $t_*\in(3.10,3.12)$ with $h_{G_1}(t_*)=h_{G_0}(t_*)$.

**Proof and reproducible certificate.** The accompanying `no83_stage29_cert.cpp` performs 160-bit MPFR directed-rounding interval arithmetic for the *discrete explicit-midpoint iterates* with step $s=10^{-4}$, after 31,000 and 31,200 steps. It encloses the discrete differences by

$$
\begin{aligned}
D_{\rm mid}(3.10)&\in[0.00013032196154450258,\;0.00013032196154450260],\\
D_{\rm mid}(3.12)&\in[-0.00026697681328268614,\;-0.00026697681328268608].
\end{aligned}
$$

These alone are enclosures of iterates rather than the exact ODE solutions. For either graph in the infinity norm, valid uniform bounds are $B=4/5$, $L=8/5$, $H=16/5$ for $\|F\|$, $\|DF\|$, and $\|D^2F\|$, respectively. For the explicit midpoint map, the local defect is at most $Cs^3$, where

$$
C=\frac{HB^2+L^2B}{6}+\frac{HB^2}{8}=0.938666666\ldots,
$$

and its Lipschitz factor is at most $q=1+Ls+L^2s^2/2$. Thus the error per graph at $n$ steps is at most $Cs^3(q^n-1)/(q-1)$. The two-graph difference has the convenient conservative error bound

$$
E(T)=2\frac CLs^2\operatorname{expm1}\!\left((L+L^2s/2)T\right).
$$

An independent audit recalculated $E(3.10)\approx1.66203\times10^{-6}<1.663\times10^{-6}$ and $E(3.12)\approx1.71646\times10^{-6}<1.717\times10^{-6}$. Hence, conservatively,

$$
\begin{aligned}
D(3.10)&\in(0.00012865,\;0.00013199),\\
D(3.12)&\in(-0.00026870,\;-0.00026526).
\end{aligned}
$$

These intervals have opposite signs. Continuity of the ODE solutions gives the claimed $t_*$. $\square$

The certificate code was reviewed line by line by an independent AI auditor; the auditor also reproduced its discrete-interval outputs using `mpmath` interval arithmetic and separately checked the signs with DOP853. The auditor did **not** run the author's MPFR executable. The Windows native execution transcript was supplied by the author, while the original native log file and its original-file hash were not supplied.

## 4. Statistical implication and limits

At the collision time $t_*$, under the observation model $Y_k=h_G(t_*)+Z_k$ with independent $Z_k\sim\mathcal N(0,\sigma^2)$, the joint observation distributions for $G_0$ and $G_1$ are identical. For the two *specified* hypotheses and equal priors, every decision rule has Bayes error at least $1/2$ (and the minimum equals $1/2$). This conclusion is limited to repeated observations at this **single** time; it says nothing about other observation times, multiple-time measurements, different resets, or general nonidentifiability.

## 5. Relation to Research No. 28 and literature

The reset and second-derivative identity were previously proposed in Research No. 28. Here the identity is proved directly in Lemma 1, independently of access to the original Appendix A; historical priority and the precise wording of that appendix have not been checked. A proposed $O(\sigma^2N^4(\log N)^8\log(2/\delta))$ sample upper bound is excluded from the theorems because it remains a proof sketch. The proposed contribution is the juxtaposition of graph-uniform guaranteed early separation with a concrete later observation collision. Literature novelty remains **INCONCLUSIVE**.

## 6. Audit, reproducibility, and publication status

- Independent AI audit dated 24 November 2026: **PASS** for Theorems 1 and 2, the analytic error argument, and source-code inspection; **not human peer review**.
- Windows author-side C++ compilation and execution on Strawberry MinGW were reported **PASS**; the chat-transcribed output is retained with its provenance limitation. Original `NO83_WINDOWS_RUN_LOG.txt` and native toolchain version output as files were **not** received; Windows provenance therefore remains **INCONCLUSIVE** at original-file level.
- References: **PARTIAL**. Novelty: **INCONCLUSIVE**. No theorem-level exhaustive literature search has been completed.
- This is the public preprint v1.0, released with the author’s authorization on 24 November 2026. **All rights reserved.**

## References (bibliographic checks remain partial)

1. R. Yonezu, *Restricted Hub-Star Phase Reconstruction*, Research No. 28, Zenodo DOI: 10.5281/zenodo.23443118; original technical Appendix A not directly compared.
2. Shandilya and Timme, *Inferring network topology from complex dynamics*, *New Journal of Physics* 13 (2011) 013004, DOI: 10.1088/1367-2630/13/1/013004 (details not independently reconciled in full).
3. *Fundamental limitations of network reconstruction from temporal data* (2017), https://pmc.ncbi.nlm.nih.gov/articles/PMC5332581/ (title/metadata not fully verified).
4. *Identifying Network Structure of Nonlinear Dynamical Systems: Contraction and Kuramoto Oscillators* (2025), arXiv:2509.13505 (precise theorem overlap not excluded).
