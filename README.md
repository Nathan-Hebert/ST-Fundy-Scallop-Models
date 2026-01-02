# ST-Fundy-Scallop-Models
Code used in Hebert et al. 2025, "Interannual variability in the length–weight relationship can disrupt the abundance–biomass correlation of sea scallop (<i>Placopecten magellanicus</i>)". ICES Journal of Marine Science.

Recommended order to run the `.R` files:
- `Survey Data Preprocessing.R`
- `Prediction Data Preprocessing.R`
- `Modelling + Predictions.R`
- `Second Sensitivity Analysis.R`

`Modelling + Predictions.R` and `Second Sensitivity Analysis.R` depend on the processed survey data from `Survey Data Preprocessing.R` and the processed environmental data from `Prediction Data Preprocessing.R`.

Survey data used in this study are available from the Government of Canada Open Data Portal ([shell height frequency data](https://open.canada.ca/data/en/dataset/ecc09d98-56ed-4a27-ad62-5c3714a1d9b4), [joint meat weight–shell height data](https://open.canada.ca/data/en/dataset/65d32794-2d81-4682-b0ea-8d8bbe907a58)).
