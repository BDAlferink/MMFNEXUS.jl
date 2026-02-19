```@meta
CurrentModule = MMFNEXUS
```

# Theory
MMF NEXUS is an algorithm for the identification of cosmic web environments. Clusters (nodes), filaments, sheets (walls), and voids are found in scale space to preserve the multiscale nature of the large scale distributions of matter. This packages concerns specifically the NEXUS+ procedure which computes the filtering in logarithmic space. 

The algorithm is outlined by the following steps:
1) Apply a gaussian filter with scale $R_n$ to the input field.
2) Compute the eigenvalues of the Hessian matrix of the filtered field.
3) Compute signatures for each environment for each voxel of the density field.
4) Repeat the steps above for a number of different filtering scales.
5) Find the optimal environment signature, independant of scale.
6) Use a detection threshhold based on physical criteria.

The environments are found hyrarchically, in the sense that first nodes are detected, then filaments, then walls, and anything not identified is considered a void. This is done as, for example, typically nodes have both a strong node and filament siganture, etc.

Below we elaborate on each step. This theory section is paraphrasiing the detailed explanation found in the original [NEXUS publication](https://academic.oup.com/mnras/article/429/2/1286/1038906).

### Filtering
The gaussian filter of width $R_n$ is applied in Fourier space. First for node identification, this is done on the density field $f(\bm{x})$,
```math
f_n(\bm{x}) = \int \frac{d^3k}{(2\pi)^3}e^{-k^2R_n^2/2} \hat{f}(\bm{k})e^{i\bm{k}\cdot\bm{x}},
```
where $\hat{f}(\bm{k})$ is the fourier transform of the input field and $f_n(\bm{x})$ is the filtered field. This is used for the input field when calculating the node signatures. They turn out to perform much worse when using the NEXUS+ procedure of filtering in logarithmic space.

For the signature calculations of the other environments, the log of the field f is first computed,
```math
g = \log_{10}{f},
```
and then smoothed with a gaussian filter,
```math
g_n(\bm{x}) = \int \frac{d^3k}{(2/pi)^3}e^{-k^2R_n^2/2} \hat{g}(\bm{k})e^{i\bm{k}\cdot\bm{x}},
```
where $\hat{g}(\bm{k})$ is the fourier transform of the logarithm of the input field and $g_n(\bm{x})$ is the logarithmic filtered field.
To obtain the lograithmically filtered field, we exponantiate and renormalize the resulting field such that the mean remains the same as the input field,
```math
f_n(\bm{x}) = C_n 10^{g_n},
```
where $C_n$ is the renormalization constant.

Usually, filtering is done over a range of scales naturally seperated as:
```math
R_n = (\sqrt{2})^n R_0,
```
where $R_0$ is the minimum filtering scale and n is choosen to be a range of integers. Typical values are $R_0 = 0.5 Mpc/h$ and $n = 0, \ldots, 6$.

### Signature calculations
With this we can compute the Hessian of the the field for each voxel. The hessian is defined as,
```math
\begin{equation}
    H_{ij}(\bm{x},R_n) = R_n^2 \frac{\partial^2 f_n(\bm{x})}{\partial x_i \partial x_j}.
\end{equation}
```
This is also computed in Fourier space.
the Hessian eigenvalues are calculated as,
```math
\det(H(\bm{x},R_n) - \lambda I) = 0,
```
where $\lambda$ are the eigenvalues.
The eigenvalues are ordered as $\lambda_1 \leq \lambda_2 \leq \lambda_3$.

With this, we find the morphological signature for each environment,
```math
\mathcal{S}(\bm{x}, R_n) = \begin{cases}
\frac{\lambda_3^2}{|\lambda_1|} \theta(-\lambda_1) \theta(-\lambda_2) \theta(-\lambda_3) & \text{node}, \\
\frac{\lambda_2^2}{|\lambda_1|} \Theta\left(1 - \frac{|\lambda_3|}{|\lambda_1|} \right) \theta(-\lambda_1) \theta(-\lambda_2) & \text{filament}, \\
|\lambda_1| \Theta\left(1 - \frac{|\lambda_2|}{|\lambda_1|} \right) \Theta\left(1 - \frac{|\lambda_3|}{|\lambda_1|} \right) \theta(-\lambda_1) & \text{wall},
\end{cases}
```
where $\Theta(x) = x \theta(x)$ and $\theta(x)$ is the Heaviside step function.
This classification is based on the insights of [Zel'dovich 1970](https://articles.adsabs.harvard.edu/pdf/1970A%26A.....5...84Z) and [Hahn et al. 2007](https://academic.oup.com/mnras/article/375/2/489/1441686) on the geometry of the matter distribution based on the eigenvalues of the Hessian.

Then, for each point, the signature is choosen for the scale where the signature is maximum,
```math
\mathcal{S}(\bm{x}) = \max_{n} \mathcal{S}(\bm{x}, R_n),
```
which is a scale independant signature.

### Detection thresholds

For nodes, the threshold is determined by considering the fraction of valid nodes given the average density contrast contained in the identified nodes as a function of the signature threshold. This is done by identifying all clumps of adjacent voxels above a certain signature threshold as individual nodes. Then, the average density of each clump is calculated and compared to a given average density $\Delta$ often taken as $\Delta = 370$. Increasing the cluster signature threshold will monotonically increase the fraction of valid nodes from near 0 to 1. The optimal threshold is then chosen as the signature corresponding to 50% valid nodes. In addition to this, we only consider nodes that contain a mass of at least $10^{13} M_\odot$ to avoid spurious detections.

Filament and wall detection thresholds are determined by considering the change in mass contained in filaments/walls as a function of signature threshold. This is given by,
```math
\Delta M^2 = \left| \frac{dM^2}{d\log S} \right|, 
```
where $M$ is the mass contained in filaments/walls and $S$ the signature threshold. The optimal signature threshold is then chosen as the signature corresponding to the peak in $\Delta M^2$. An additional requirement for the minimum volume of filaments and walls is chosen to be 1 $(\text{Mpc}/h)^3$ and 2 $(\text{Mpc}/h)^3$, respectively.

### Variations
There are a number of ``other`` routines that find signatures of different fields such the tidal field, the velocity divergence, the velocity shear, the logdensity field (this is not the same as NEXUS+ which only applies the filter in log space). As of now, none of those are implemented in this package but may be added in the future. More information is found in the original [NEXUS publication](https://academic.oup.com/mnras/article/429/2/1286/1038906).