
## TODO
* Put cleaned-up code required to run analyese into the respeicte folders (sorted by priority)
  * [ ] `association_testing` [Brian]
  * [x] `comparison_methods`
     * [x]  explain what annoation data frame columns are needed
     * [x] `monti` [Eva]
     * [x] `staar` [Eva]
  * [x] `phenotype_prediction` [Eva]
  * [x] `simulation` [Eva]
  * [x] `feature_importance` [Hakime]
  * [ ] `repeat analysis` (no folder yet, to be decided) [Brian] 
  * [ ] `correlation analysis` (no folder yet, to be decided)[Brian]
* Put scripts/markdowns to generate final plots int `paper_figures`
   * name according to the figure they belong to in the main text


## Cleanup

For each file:
* Remove orphaned functions/rules
  1. Remove rules not needed in pipelines
  1. Remove click commands not used in any pipeline
  1. Remove functions not used in any Click command
* Remove unused imports, global variables, command/function arguments
* Remove all comments that are not explanations of code
* Remove `resources` from snakefiles
* Convert `genopheno` imports to corresponding `deeprvat` imports
* Reformat files with black/snakefmt. Use formatR::tidy_source/formatR::tidy_dir for r code
