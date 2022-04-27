// Copyright 2014 University of Notre Dame
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

using Landis.Utilities;

namespace Landis.Library.Succession.DemographicSeeding
{
    /// <summary>
    /// Species parameters for demographic seeding that the user can specify.
    /// </summary>
    public class SpeciesParameters
    {
        private double minSeedsProduced;
        private double maxSeedsProduced;
        private double seedMass;
        private double slwMax;
        private double fracFol;
        private double frActWd;
        private double seedCalibration;
        public double[] DispersalParameters { get; private set; }
        public double[] EmergenceProbabilities { get; private set; }
        public double[] SurvivalProbabilities { get; private set; }
        public double[] MaxSeedBiomass { get; private set; }

        // Indexes for dispersal parameters
        public const int DoubleExponential_Mean1   = 0;
        public const int DoubleExponential_Mean2   = 1;
        public const int DoubleExponential_Weight1 = 2;

        //---------------------------------------------------------------------

        public SpeciesParameters()
        {
            DispersalParameters = new double[3];
            EmergenceProbabilities = new double[Model.Core.Ecoregions.Count];
            SurvivalProbabilities = new double[Model.Core.Ecoregions.Count];
            MaxSeedBiomass = new double[Model.Core.Ecoregions.Count];
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Minimum number of seeds produced by an occupied cell in a year.
        /// </summary>
        public double MinSeedsProduced
        {
            get
            {
                return minSeedsProduced;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "Minimum seeds produced must be = or > 0");
                minSeedsProduced = value;
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Maximum number of seeds produced by an occupied cell in a year.
        /// </summary>
        public double MaxSeedsProduced
        {
            get
            {
                return maxSeedsProduced;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "Maximum seeds produced must be = or > 0");
                maxSeedsProduced = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Average seed mass
        /// </summary>
        public double SeedMass
        {
            get
            {
                return seedMass;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "Seed mass must be = or > 0");
                seedMass = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Specific leaf weight max
        /// </summary>
        public double SLWmax
        {
            get
            {
                return slwMax;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "SLWmax must be = or > 0");
                slwMax = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Fraction Active Wood
        /// </summary>
        public double FrActWd
        {
            get
            {
                return frActWd;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "FrActWd must be = or > 0");
                frActWd = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Fraction Foliage
        /// </summary>
        public double FracFol
        {
            get
            {
                return fracFol;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "FracFol must be = or > 0");
                fracFol = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Seed Calibration factor
        /// </summary>
        public double SeedCalibration
        {
            get
            {
                return seedCalibration;
            }
            set
            {
                if (value < 0)
                    throw new InputValueException(value.ToString(),
                                                  "SeedCalibration must be = or > 0");
                seedCalibration = value;
            }
        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Mean1 parameter for 2-component exponential dispersal kernel
        /// </summary>
        public double DispersalMean1
        {
            get
            {
                return DispersalParameters[DoubleExponential_Mean1];
            }
            set
            {
                if (value < 0.0)
                    throw new InputValueException(value.ToString(),
                                                  "Dispersal Mean1 must be = or > 0.0");
                DispersalParameters[DoubleExponential_Mean1] = value;
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Mean2 parameter for 2-component exponential dispersal kernel
        /// </summary>
        public double DispersalMean2
        {
            get
            {
                return DispersalParameters[DoubleExponential_Mean2];
            }
            set
            {
                if (value < 0.0)
                    throw new InputValueException(value.ToString(),
                                                  "Dispersal Mean2 must be = or > 0.0");
                DispersalParameters[DoubleExponential_Mean2] = value;
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Weight1 parameter for 2-component exponential dispersal kernel
        /// </summary>
        public double DispersalWeight1
        {
            get
            {
                return DispersalParameters[DoubleExponential_Weight1];
            }
            set
            {
                if (value < 0.0)
                    throw new InputValueException(value.ToString(),
                                                  "Dispersal Weight1 must be = or > 0.0");
                DispersalParameters[DoubleExponential_Weight1] = value;
            }
        }
    }
}
