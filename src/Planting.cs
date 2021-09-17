using Landis.Core;
using System.Collections;
using System.Collections.Generic;
using Landis.SpatialModeling;

namespace Landis.Library.Succession
{
    /// <summary>
    /// The form of reproduction that represents the planting of trees after
    /// a harvest.
    /// </summary>
    public class Planting
        : FormOfReproduction, IPlanting
    {
        /// <summary>
        /// A list of species to be planted
        /// </summary>
        public class SpeciesList
        {
            private BitArray bitArray;
            
            // JSF - Convert dictionary to float
            private Dictionary<ISpecies, uint> plantDict;

            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            public SpeciesList(IEnumerable<ISpecies> speciesList,
                               ISpeciesDataset speciesDataset)
            {
                bitArray = new BitArray(speciesDataset.Count);
                if (speciesList != null) {
                    foreach (ISpecies species in speciesList) {
                        bitArray.Set(species.Index, true);
                    }
                }
            }

            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            public SpeciesList(IEnumerable<ISpecies> speciesList,
                               ISpeciesDataset speciesDataset,
                               Dictionary<ISpecies, uint> plantNum)
            {
                bitArray = new BitArray(speciesDataset.Count);
                plantDict = new Dictionary<ISpecies, uint>();
                if (speciesList != null)
                {
                    foreach (ISpecies species in speciesList)
                    {
                        bitArray.Set(species.Index, true);
                        plantDict[species] = plantNum[species];
                    }

                }
            }

            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            internal BitArray AsBitArray
            {
                get {
                    return bitArray;
                }
            }

            internal Dictionary<ISpecies, uint> AsPlantDict
            {
                get
                {
                    return plantDict;
                }
            }
        }

        //---------------------------------------------------------------------

        public Planting()
            : base()
        {
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Are the conditions necessary for planting a species at site
        /// satified?
        /// </summary>
        protected override bool PreconditionsSatisfied(ISpecies   species,
                                                       ActiveSite site)
        {
            //return true; 
            return Reproduction.PlantingEstablish(species, site);
        }

        //---------------------------------------------------------------------

        void IPlanting.Schedule(SpeciesList speciesList,
                                ActiveSite  site)
        {
            SelectedSpecies[site].Or(speciesList.AsBitArray);
            if (speciesList.AsPlantDict != null)
            {
                foreach(var plant in speciesList.AsPlantDict)
                {
                    PlantingList[site].Add(plant.Key, plant.Value);
                }
            }
        }
    }
}
