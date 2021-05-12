using System;
using System.Diagnostics;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Troschuetz.Random;
using Landis.Core;
using Landis.SpatialModeling;
using Accord.Statistics.Distributions.Univariate;
using MathNet.Numerics;
using MathNet.Numerics.Integration;


namespace Landis.Library.Succession.DensitySeeding
{
    public enum Dispersal_Model { DOUBLE_EXPONENTIAL, TWODT };
    public enum Seed_Model { FIXED, UNIFORM, BIOMASS, DENSITY };         // a enumerable that identifies a particular seed production model 
    public enum Dispersal_Type { STATIC, DYNAMIC };    // either “STATIC”, meaning that the dispersal model does not depend on circumstances that can change during a simulation, or “DYNAMIC”

    public struct Species
    {
        public int index;
        public double max_seed;
        public double min_seed;
        public double seed_mass;
        public double SLWmax;
        public double seedCalibration;
        public double[] dispersal_parameters;
        public int shade_tolerance;
        // public double leaf_area;
        public int reproductive_age;
        public int reproductive_age_steps;
        public double max_dispersal_distance;
        public int max_dispersal_distance_pixels;
    };

    public class DensitySeedMap
    {
        public int map_width_pixels;                       // map width (in pixels)
        public int map_height_pixels;                      // map height (in pixels)
        public int num_species;                            // total number of species in the simulation
        public int num_time_steps;                         // the number of succession time steps to loop over
        public int num_ecoregions;                         // total number of ecoregions in the simulation
        public int max_age_steps;                          // maximum age allowed for any species, in years
        public int time_step;                              // the length of a succession time step (in years)
        public double pixel_size;                          // size of a raster pixel (in meters)
        public int mc_draws;                               // the number of Monte Carlo draws to use when estimating dispersal probabilities
        public double max_leaf_area;                       // maximum projected seedling canopy area that can be supported in a cell
        public int cohort_threshold;                       // minimum number of trees to qualify as a cohort under the default succession model
        public double seedling_leaf_area;                  // Area (m2) that each seedling occupies (all species)
        public double min_cohort_prop;                      // Minimum Proportion of max biomass that counts as a cohort

        public int rand_seed;                              // seed for random number generator

        public int[][] ecoregion;                          // a 2 dimensional integer array (raster map) that includes the entire study area and defines an ecoregion for each cell 
        public int[][] shade;                              // a 2 dimensional integer array (raster map) that describes the understory shade available in each map cell.  May take values from 0 to 5.  Has the same dimensions as ecoregion.;
        public bool[][][] reproductive;                     // a 3 dimensional integer array in which the s-th slice along the first dimension is a map that describes the presence (= 1) or absence (= 0) of a reproductive-aged population on each cell for species s.  Dimensions are (n_spp, mapdim[1], mapdim[2]).  Represented as an integer instead of Boolean to preserve option to change to a count in the future
        public int[][][][] cohorts;                        // a 4 dimensional integer array containing a map of each species-age cohort

        public double[][][][] seed_production;                  // a 4 dimensional double float point array in which the s-th slice along the first dimension is a map that describes the number of seeds from species s that are produced in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).  
        public double[][][][] seed_shadow;                     // a 4 dimensional integer array in which the s-th slice along the first dimension is a map that describes the number of seeds from species s that arrived in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).  
        public double[][][][] n;
        public double[][][][] v;
        public double[][][][] seed_emergence;
        public int[][][] seedlings;                         // a 3 dimensional integer array in which the s-th slice along the first dimension is a map that describes the number of newly established trees from species s in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).  
        public double[][][] fol_mass;                        // a 3 dimensional double float point array in which the s-th slice along the first dimension is a map that describes the mature foliage mass for species s in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).
        public double[][][] density_seeds;                        // a 3 dimensional double float point array in which the s-th slice along the first dimension is a map that describes the number of seeds produced for species s in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).
        public double[][] max_seed_biomass;

        public double[][] survival_probability;
        private Dictionary<double, double>[] dispersal_probability;

        public Dispersal_Model dispersal_model;
        public Seed_Model seed_model;
        public Dispersal_Type dispersal_type;

        public Species[] all_species;

        //----Density seed dispersal----
        public int[][][] seedProduction;
        public int[][][] seedDispersal;
        public double[][][] emergence_probability;          // a 3 dimensional double float point array in which the s-th slice along the first dimension is a map that describes the probabiltiy of establishment for species s in each map cell during the current time step. Dimensions are (n_spp, mapdim[1], mapdim[2]).

        float alpha;
        int prob_n;

        public int[] numcell_all_list = new int[Model.Core.Species.Count];
        public int[] hilmt4spec;
        public int[] wghts;
        public double[][][] prob4square;
        public double[][] cdf4species;
        double[] minprob4spe = new double[Model.Core.Species.Count];
        double constsca;

        public DensitySeedMap(int timeStep)
        {
            int maxCohortAge;  // maximum age allowed for any species, in years

            map_width_pixels = Model.Core.Landscape.Columns;
            map_height_pixels = Model.Core.Landscape.Rows;
            num_species = Model.Core.Species.Count;
            num_time_steps = (Model.Core.EndTime - Model.Core.StartTime) / timeStep;
            time_step = timeStep;
            num_ecoregions = Model.Core.Ecoregions.Count;

            maxCohortAge = 0;
            foreach (ISpecies species in Model.Core.Species)
                if (species.Longevity > maxCohortAge)
                    maxCohortAge = species.Longevity;

            max_age_steps = maxCohortAge / timeStep;

            // Initialize seed production and seed dispersal arrays
            hilmt4spec = new int[Model.Core.Species.Count];
            prob4square = new double[Model.Core.Species.Count][][];
            cdf4species = new double[Model.Core.Species.Count][];
            seedProduction = new int[num_species][][];
            seedDispersal = new int[num_species][][];
            emergence_probability = new double[num_species][][];

            for (int s = 0; s < num_species; s++)
            {
                seedProduction[s] = new int[map_width_pixels][];
                seedDispersal[s] = new int[map_width_pixels][];
                emergence_probability[s] = new double[map_width_pixels][];
                for (int x = 0; x < map_width_pixels; x++)
                {
                    seedProduction[s][x] = new int[map_height_pixels];
                    seedDispersal[s][x] = new int[map_height_pixels];
                    emergence_probability[s][x] = new double[map_height_pixels];
                }
            }

        }



        public int initialize(float speciesMaxDist, float speciesEffectiveDist, int speciesID)
        {
            Debug.Assert(speciesMaxDist > 0);


            double[] alpha_x = { 9.233413476451586, 11.228872242412663, 13.062240779188071, 14.794149222537211, 16.454745203680105, 18.061636840199071, 19.626177395384239, 21.156198165839982, 22.657373309062930 };

            for (int i = 0; i < 9; i++)
            {
                int n = i + 2;

                prob_n = n;

                alpha = (float)(alpha_x[i] / Math.Pow(speciesMaxDist, 1.0 / n));

                double x_val = (alpha * Math.Pow(speciesEffectiveDist, 1.0 / n));

                int ret_id = check_prob(x_val, n);

                if (ret_id == 0)
                    break;

                if (i == 8 && ret_id == -1)
                {
                    //FIXME JSF Error message
                    //printf("\n\nplease check effective distance and max distance of the %d th species in SpeciesAtrributes.dat.\n");
                    //printf("something might be wrong.\n within the effective distance, the cumulated probability is not greater than 0.95\n");
                    //Sleep(10000);
                }
            }

            int tmp_offset = (int)Math.Round(speciesMaxDist / Model.Core.CellLength);
            int offset = tmp_offset >= 1 ? tmp_offset : 1;
            int numcellside = 2 * offset + 1;
            int numcell_all = numcellside * numcellside;

            var binom_tro = new Troschuetz.Random.Distributions.Discrete.BinomialDistribution(0.17, 2100);
            int trobinom = binom_tro.Next();

            basefunction basefunc;
            basefunc = new basefunction(alpha, prob_n, wghts);

            // Test
            Random rand = new Random();
            var binomial_distribution = new BinomialDistribution(2100, 0.17);
            int testbinom = binomial_distribution.Generate(rand);

            numcell_all_list[speciesID] = numcell_all;

            hilmt4spec[speciesID] = offset;

            prob4square[speciesID] = new double[numcellside][];
            //rob4square[speciesID][0] = new double [numcell_all];

            for (int i = 0; i < numcellside; i++)
            {

                prob4square[speciesID][i] = new double[numcell_all];

            }

            cdf4species[speciesID] = new double[numcell_all];

            return numcellside;
        }

        public void DispereSiteSeeds(ISpecies species, ActiveSite site)
        {
            // clock_t start_s=clock();
            Random rand = new Random();

            int x = site.Location.Column - 1;
            int y = site.Location.Row - 1;
            int s = species.Index;
            int offset = hilmt4spec[s];

            int totalSeeds = seedProduction[s][x][y];
            int sumofDispersal = 0;

            if (totalSeeds > 0)
            {

                for (int i = -offset; i <= offset; i++)
                {
                    int spread_row = y + i;

                    if (spread_row < 0 || spread_row >= map_width_pixels)
                        continue;

                    int prob4square_row_id = i + offset;

                    for (int j = -offset; j <= offset; j++)
                    {
                        int spread_col = x + j;

                        if (spread_col < 0 || spread_col >= map_height_pixels)
                            continue;

                        int prob4square_col_id = j + offset;

                        var binomial_distribution = new BinomialDistribution(totalSeeds, prob4square[s][prob4square_row_id][prob4square_col_id]);
                        int number = (int)binomial_distribution.Generate(rand);
                        seedDispersal[s][spread_row][spread_col] += number;
                        sumofDispersal += number;
                    }
                }
            }
            //seedDispersal[s][x][y] = Math.Max((totalSeeds - sumofDispersal), 0);

        }

        //---------------------------------------------------------------------
        public void CheckEstablishment(ISpecies species, ActiveSite site)
        {
            Random rand = new Random(); 
            
            int x = site.Location.Column - 1;
            int y = site.Location.Row - 1;
            int s = species.Index;

            int totalSeeds = seedDispersal[s][x][y];
            double establishProb = emergence_probability[s][x][y];

            var binomial_distribution = new BinomialDistribution(totalSeeds, establishProb);
            seedDispersal[s][x][y] = binomial_distribution.Generate(rand);
        }

        //---------------------------------------------------------------------
        public int check_prob(double x_val, int n)
        {

            int demon = 1;

            int[] weight = new int[n];
            double[] x_list = new double[n];

            for (int i = n - 1; i >= 1; i--)
                demon *= i;

            constsca = 1 / (2 * Math.PI * demon);

            weight[0] = 1;
            x_list[n - 1] = 1;

            for (int i = n - 1; i >= 1; i--)
            {
                weight[n - i] = weight[n - i - 1] * i;
                x_list[i - 1] = x_list[i] * x_val;
            }

            double term = 0;

            for (int i = 0; i < n; i++)
                term += x_list[i] * weight[i];

            double prob = 1 - Math.Exp(-x_val) * term / demon;

            if (prob > 0.94)
                return 0;
            else
                return -1;
        }

        public double calculateSeedingProbability(float x_cor1, float x_cor2, float y_cor1,
            float y_cor2, int relative_row, int relative_col)
        {
            Debug.Assert(x_cor1 < x_cor2);
            Debug.Assert(y_cor1 < y_cor2);

            float xcor1 = Math.Abs(x_cor1);
            float xcor2 = Math.Abs(x_cor2);
            float ycor1 = Math.Abs(y_cor1);
            float ycor2 = Math.Abs(y_cor2);

            double prob = 0;

            if (relative_row != 0 && relative_col != 0)
            {
                float x1 = Math.Min(xcor1, xcor2);
                float x2 = Math.Max(xcor1, xcor2);
                float y1 = Math.Min(ycor1, ycor2);
                float y2 = Math.Max(ycor1, ycor2);

                prob = cal_rect_prob(x2, y2) - cal_rect_prob(x1, y2) - cal_rect_prob(x2, y1) + cal_rect_prob(x1, y1);
            }
            else
            {
                if (relative_row == 0 && relative_col == 0)
                {
                    prob = cal_rect_prob(x_cor2, y_cor2) * 4;
                }
                else if (relative_row == 0)
                {
                    Debug.Assert(relative_col != 0);

                    float y1 = Math.Min(ycor1, ycor2);
                    float y2 = Math.Max(ycor1, ycor2);

                    prob = (cal_rect_prob(x_cor2, y2) - cal_rect_prob(x_cor2, y1)) * 2;
                }
                else
                {
                    Debug.Assert(relative_row != 0 && relative_col == 0);

                    float x1 = Math.Min(xcor1, xcor2);
                    float x2 = Math.Max(xcor1, xcor2);

                    prob = (cal_rect_prob(x2, y_cor2) - cal_rect_prob(x1, y_cor2)) * 2;
                }
            }

            Debug.Assert(prob >= -1.0e-14 && prob <= 1 + 1.0e-14);

            return prob;

        }


        public void calculateSquareProbability(int speciesID, int numCellSide)
        {
            wghts = new int[prob_n];
            wghts[0] = 1;

            for (int i = prob_n - 1; i >= 1; i--)
                wghts[prob_n - i] = wghts[prob_n - i - 1] * i;

            //================================================================
            basefunction basefunc;
            basefunc = new basefunction(alpha, prob_n, wghts);

            int numcell_all = numCellSide * numCellSide;

            int hilmt = hilmt4spec[speciesID];
            int lwlmt = -hilmt;
            float cell_size = Model.Core.CellLength;

            double sum = 0;
            int ii = 0;
            int jj = 0;
            for (int relative_col = hilmt; relative_col >= lwlmt; relative_col--)
            {
                float y_cor1 = (relative_col - 0.5f) * cell_size;
                float y_cor2 = (relative_col + 0.5f) * cell_size;

                for (int relative_row = lwlmt; relative_row <= hilmt; relative_row++)
                {
                    float x_cor1 = (relative_row - 0.5f) * cell_size;
                    float x_cor2 = (relative_row + 0.5f) * cell_size;

                    if (prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col] < 1.0e-14)
                        prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col] = calculateSeedingProbability(x_cor1, x_cor2, y_cor1, y_cor2, relative_row, relative_col);

                    if (prob4square[speciesID][relative_row - lwlmt][hilmt + relative_col] < 1.0e-14)
                        prob4square[speciesID][relative_row - lwlmt][hilmt + relative_col] = prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col];

                    if (prob4square[speciesID][hilmt - relative_col][relative_row - lwlmt] < 1.0e-14)
                        prob4square[speciesID][hilmt - relative_col][relative_row - lwlmt] = prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col];

                    if (prob4square[speciesID][hilmt + relative_col][relative_row - lwlmt] < 1.0e-14)
                        prob4square[speciesID][hilmt + relative_col][relative_row - lwlmt] = prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col];

                    sum += prob4square[speciesID][relative_row - lwlmt][hilmt - relative_col];
                    jj += 1;
                }
                ii += 1;
            }

            // FILE* fp = fopen("o.txt", "w");

            // for(int relative_col=hilmt; relative_col>=lwlmt; relative_col--)
            // {
            // 	const float y_cor1 = (relative_col - 0.5f) * cell_size;
            // 	const float y_cor2 = (relative_col + 0.5f) * cell_size;

            // 	for(int relative_row=lwlmt; relative_row<=hilmt; relative_row++)
            // 	{
            // 		const float x_cor1 = (relative_row - 0.5f) * cell_size;
            // 		const float x_cor2 = (relative_row + 0.5f) * cell_size;

            // 		fprintf(fp, "%8.3lf", prob4square[spec_id][relative_row - lwlmt][hilmt - relative_col] * 1000);
            // 	}

            // 	fprintf(fp,"\n");
            // }

            // fclose(fp);

            double minprob = prob4square[speciesID][0][0];
            double maxprob = prob4square[speciesID][hilmt][hilmt];

            minprob4spe[speciesID] = minprob;
            // printf("min = %e, max = %e sum = %lf\n", minprob, maxprob, sum);

            for (int i = 0; i < numcell_all; i++)
            {
                //Debug.Assert(minprob <= prob4square[speciesID][0][i]);
                //Debug.Assert(maxprob >= prob4square[speciesID][0][i]);
            }



        }

        public double cal_rect_prob(float x, float y)
        {
            Debug.Assert(x > 0 && y > 0);

            // Test
            Random rand = new Random();
            var binomial_distribution = new BinomialDistribution(2100, 0.17);
            int testbinom = binomial_distribution.Generate(rand);

            double prob_b = prob_base(x, y);

            double prob = 0.25 - constsca * prob_b;

            return prob;
        }

        public double prob_base(float x, float y)
        {
            basefunction basefunc;
            basefunc = new basefunction(alpha, prob_n, wghts);

            basefunc.set_onecor(x);

            double theta1 = Math.Atan(y / x);
            double integrate1 = DoubleExponentialTransformation.Integrate(basefunc.operate, 0.0, theta1, 1e-8);


            basefunc.set_onecor(y);
            double theta2 = Math.PI / 2 - theta1;
            double integrate2 = DoubleExponentialTransformation.Integrate(basefunc.operate, 0.0, theta2, 1e-8);

            double ret_val = integrate1 + integrate2;

            return ret_val;
        }

        class basefunction
        {

            private static double alpha;
            private static double one_cor; //one coordinate: x or y
            private static int n;
            private static int[] wghts;
            private static double[] z_list;
            public basefunction(double alpha_in, int prob_n, int[] wghts_in)
            {
                alpha = alpha_in;
                n = prob_n;
                wghts = wghts_in;

                z_list = new double[n];
            }


            public void set_onecor(double one_cor_in)
            {
                one_cor = one_cor_in;
            }


            public Func<double, double> operate = delegate (double theta)
            {

                double z_val = Math.Pow(one_cor / Math.Cos(theta), 1.0 / n) * alpha;

                z_list[n - 1] = 1;

                for (int i = n - 1; i >= 1; i--)
                    z_list[i - 1] = z_list[i] * z_val;

                double term = 0;
                for (int i = 0; i < n; i++)
                    term += z_list[i] * wghts[i];

                double retval = Math.Exp(-z_val) * term;

                return retval;
            };
        }

        private bool isInside(int x, int y) // check if the coordinates are inside the map
        {
            return (x >= 0 && y >= 0 && x < map_width_pixels && y < map_height_pixels);
        }

    };
}