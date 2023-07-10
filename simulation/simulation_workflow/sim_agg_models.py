import torch
import torch.nn as nn
from scipy.stats import beta
import numpy as np


# TODO: Make weighting functino more flexible


class GenoMAFWeight(nn.Module):
    def __init__(self, maf_dim: int = 0):
        # , sigmoid: bool = False):
        super().__init__()

        self.maf_dim = maf_dim

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        print(f"X.shape: {X.shape}")
        t = torch.abs(
            torch.log10(X[:, :, self.maf_dim, :])
        )  # select dim 0 assuming that this is the dimension with MAF annotation
        t[t == float("Inf")] = 0

        burden_per_gene = torch.sum(t, dim=2)  # gene_scores
        print(f"burden shape: {burden_per_gene.shape}")
        return burden_per_gene


class GenoMAFBeta(nn.Module):
    def __init__(self, beta_a: int = 1, beta_b: int = 25, maf_dim: int = 0):
        # , sigmoid: bool = False):
        super().__init__()
        self.a = beta_a
        self.b = beta_b
        self.maf_dim = maf_dim

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        t = beta.pdf(
            X[:, :, self.maf_dim, :], self.a, self.b
        )  # select dim 0 assuming that this is the dimension with MAF annotation

        # import ipdb; ipdb.set_trace()

        burden_per_gene = torch.sum(torch.Tensor(t), dim=2)  # gene_scores

        return burden_per_gene


#### ADDING ANNOTATION INCLUDED FUNCTION HERE || HAKIME


class GenoAnnotWeight(nn.Module):
    def __init__(
        self,
        weights: np.array,
        maf_dim: int = 0,
        anno_agg_method: str = "sum",
        var_agg_method: str = "sum",
        first_agg_dim: str = "variants",
    ):
        # , sigmoid: bool = False):
        super().__init__()
        ## pick aggregation function here

        # input tensor dims : samples x genes x annotations x variants

        self.maf_dim = maf_dim
        self.anno_agg_method = anno_agg_method
        self.var_agg_method = var_agg_method

        if first_agg_dim == "annotations":
            # first aggregate annotations for each variant
            # --> vector (n_variants x 1) for each gene
            self.first_agg_dim = 2
        elif first_agg_dim == "variants":
            # first aggregate variants for each annotation
            # --> vector (1 x n_annotations) for each gene
            self.first_agg_dim = 3
        else:
            ValueError("first_agg_dim must either be 'annotations' or 'variants'")

        self.weights = torch.tensor(weights)
        if len(self.weights.shape) == 1:
            self.weights = self.weights.unsqueeze(dim=1)
        assert self.weights.shape[1] == 1

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        # maf = torch.abs(torch.log10(X[:, :, self.maf_dim, :]))
        # -log10 transform MAF
        X[:, :, self.maf_dim, :] = torch.abs(torch.log10(X[:, :, self.maf_dim, :]))
        X[:, :, self.maf_dim, :][X[:, :, self.maf_dim, :] == float("Inf")] = 0

        # TODO: check if this is really across the right dimension
        t = torch.mul(X, self.weights)

        X = aggregate_tenosr(X, self.var_agg_method, self.first_agg_dim)
        burdens_per_gene = aggregate_tenosr(X, self.anno_agg_method, 2)

        # print(f'burden shape: {burdens_per_gene.shape}')

        return burdens_per_gene


def aggregate_tenosr(t, agg_method, dim):
    top_k = 1 if t.shape[dim] == 1 else 2

    burden_choice = {
        "sum": torch.sum(t, dim=dim),
        "max": torch.max(t, dim=dim),
        "min": torch.min(t, dim=dim),
        # TODO: check if mean2 is correct implemented
        "mean2": torch.mean(torch.topk(t, top_k, dim=dim).values, dim=dim),
    }
    return burden_choice[agg_method]
