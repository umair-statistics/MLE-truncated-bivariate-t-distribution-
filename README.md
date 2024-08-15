# Maximum Likelihood Estimation of Mixing Angle θ from Truncated Bivariate t-Distribution

This repository contains R codes for obtaining the maximum likelihood estimation (MLE) of the mixing angle θ from the truncated bivariate t-distribution as proposed by [Shaw et al. (2008)](https://doi.org/10.1016/j.jmva.2007.08.006). The methods implemented here follow the procedures outlined in the paper [Umair et al. (2024)](https://doi.org/10.1177/09622802241267808) and are designed to handle truncation in bivariate distributions.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)
- [Contributing](#contributing)
- [License](#license)

## Introduction

In most pre- and post-measurement intervention studies, samples are often drawn based on cutoff points or thresholds. This selective sampling process leads to truncation in the distribution of the respective test statistics. As a result, the standard estimation techniques for distribution parameters become inappropriate because they do not account for the truncated nature of the data. It is therefore crucial to estimate the parameters of the complete distribution based on truncated samples.

[Umair et al.](https://doi.org/10.1177/09622802241267808) proposed a method for maximum likelihood estimation (MLE) of the mixing angle θ in the context of truncated bivariate t-distributions numerically. The methodology allows for accurate parameter estimation despite the truncation. This repository provides an R implementation of the proposed method, enabling researchers to perform MLE for the mixing angle θ in such scenarios.

## Installation

To use the R scripts in this repository, you can clone the repository to your local machine:

```bash
git clone https://github.com/umair-statistics/MLE-truncated-bivariate-t-distribution-.git
```
Once cloned, navigate to the directory where the MLEs.R file is located:

```bash
cd MLE-truncated-bivariate-t-distribution-
```

You will need to have R installed on your system along with some required packages. These can be installed via CRAN:

```R
devtools::install_github('umair-statistics/RTM')
install.packages(reshape2)
install.packages(ggplot2)
install.packages(ggpubr)
install.packages(latex2exp)
```

## Usage

Once you have cloned the repository and installed the required packages, you can run the scripts directly in R.

1. **Data Preparation:** Ensure your data is appropriately truncated and scaled to follow bivariate t-distribution.
2. **Running the Estimation:** Use the provided functions to estimate the mixing angle θ.

For a quick start, refer to the `MLEs.R` script in the repository.

## References

Shaw, W. T., Lee, W. M., & Wong, K. F. (2008). Maximum likelihood estimation of the mixing angle in a truncated bivariate t-distribution. *Journal of Multivariate Analysis*, 99(6), 1276-1287. DOI: [10.1016/j.jmva.2007.08.006](https://doi.org/10.1016/j.jmva.2007.08.006)

Umair, M., Khan, M., & Olivier, J. (2024). Accounting for regression to the mean under the bivariate t-distribution. *Statistical Methods in Medical Research*, [09622802241267808](https://doi.org/10.1177/09622802241267808).

## Contributing

Contributions are welcome! If you have improvements or fixes, feel free to fork this repository, make your changes, and submit a pull request.

## License
[![License (MIT)](https://img.shields.io/badge/license-MIT-blue.svg?style=plastic)](http://opensource.org/licenses/MIT)

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

Feel free to modify this template to better suit your specific needs.
