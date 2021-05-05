using Landis.Core;
using Landis.SpatialModeling;
using log4net;
using Seed_Dispersal;
using System.Reflection;
using System.Collections.Generic;
using System;
using System.IO;
using Landis.Utilities;
using System.Linq;
using Landis.Utilities.Diagnostics;
using System.Diagnostics;
using MathNet.Numerics;
using MathNet.Numerics.Integration;

namespace Landis.Library.Succession.DensitySeeding
{    
    public class Algorithm
    {
        private DensitySeeding.DensitySeedMap seedingData;
        private int timeAtLastCall = -99999;

        public Algorithm(int successionTimestep)
        {
            int numTimeSteps;  // the number of succession time steps to loop over
            int maxCohortAge;  // maximum age allowed for any species, in years

            numTimeSteps = (Model.Core.EndTime - Model.Core.StartTime) / successionTimestep;
            maxCohortAge = 0;
            foreach (ISpecies species in Model.Core.Species)
                if (species.Longevity > maxCohortAge)
                    maxCohortAge = species.Longevity;

            int max_age_steps = maxCohortAge / successionTimestep;


           
            seedingData = new DensitySeeding.DensitySeedMap(successionTimestep);

            for (int s = 0; s < Model.Core.Species.Count; s++)
            {
                ISpecies species = Model.Core.Species[s];
                int numCellSide = seedingData.initialize(species.MaxSeedDist, species.EffectiveSeedDist, species.Index);
                seedingData.calculateSquareProbability(species.Index, numCellSide);

            }

        }




        //---------------------------------------------------------------------
        //---------------------------------------------------------------------

        /// <summary>
        /// Seeding algorithm: determines if a species seeds a site.
        /// <param name="species"></param>
        /// <param name="site">Site that may be seeded.</param>
        /// <returns>true if the species seeds the site.</returns>
        public void DoesSpeciesSeedSite(ISpecies species,
                                        ActiveSite site, out bool established, out double seedlingProportion)
        {
            // Is this the first site for the current timestep?
            if (Model.Core.CurrentTime != timeAtLastCall)
            {
                SimulateOneTimestep();
                //FIXME JSF
                //WriteOutputMaps();
                //WriteProbabilities();  // Right now probabilities will not vary by timestep, so no need to output multiple times.  File is written during initialization.
            }
            timeAtLastCall = Model.Core.CurrentTime;

            int x = site.Location.Column - 1;
            int y = site.Location.Row - 1;
            int s = species.Index;
            seedlingProportion = seedingData.seedDispersal[s][x][y];

            if (seedlingProportion >= 1.0)
            {
                established = true;
            }
            else
                established = false;
        }


        protected void SimulateOneTimestep()
        {
            // Reset mapped values and update seed production
            foreach (ActiveSite site in Model.Core.Landscape)
            {
                int x = site.Location.Column - 1;
                int y = site.Location.Row - 1;
                foreach (ISpecies species in Model.Core.Species)
                {
                    int s = species.Index;
                    // Initialize to zero
                    seedingData.seedDispersal[s][x][y] = 0;
                    seedingData.seedProduction[s][x][y] = 0;

                    // Update number of seeds 
                    seedingData.seedProduction[s][x][y] = Reproduction.DensitySeeds(species, site);
                    // Update the EmergenceProbability to be equal to EstablishmentProbability from Succession extension
                    seedingData.emergence_probability[s][x][y] = Reproduction.EstablishmentProbability(species, (ActiveSite)site);
                }
            }

            // Disperse seeds
            foreach (ActiveSite site in Model.Core.Landscape)
            {
                int x = site.Location.Column - 1;
                int y = site.Location.Row - 1;

                foreach (ISpecies species in Model.Core.Species)
                {
                    int s = species.Index;

                    if (Reproduction.MaturePresent(species, site))
                    {
                        // This will cause SimOneTimestep to consider the
                        // species as reproductive at this site.
                        seedingData.DispereSiteSeeds(species, site);
                    }
                }
            }

            // Calculate survival
            foreach (ActiveSite site in Model.Core.Landscape)
            {
                int x = site.Location.Column - 1;
                int y = site.Location.Row - 1;

                foreach (ISpecies species in Model.Core.Species)
                {
                    int s = species.Index;

                    if (seedingData.seedDispersal[s][x][y] > 0)
                    {
                        seedingData.CheckEstablishment(species, site);
                    }
                }
            }
        }

    }
}
