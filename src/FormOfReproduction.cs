using Landis.Core;
using System.Collections;
using System.Collections.Generic;
using Landis.SpatialModeling;

namespace Landis.Library.Succession
{
    /// <summary>
    /// A base class for forms of reproduction.
    /// </summary>
    public abstract class FormOfReproduction
        : IFormOfReproduction
    {
        //private static Species.IDataset speciesDataset;
        private static ISpeciesDataset speciesDataset;

        //---------------------------------------------------------------------

        static FormOfReproduction()
        {
            speciesDataset = Model.Core.Species;
        }

        //---------------------------------------------------------------------

        private ISiteVar<BitArray> selectedSpecies;
        private ISiteVar<Dictionary<ISpecies, uint>> plantingList;

        //---------------------------------------------------------------------

        /// <summary>
        /// The species that have been selected for this form of reproduction
        /// at each active site.
        /// </summary>
        protected ISiteVar<BitArray> SelectedSpecies
        {
            get {
                return selectedSpecies;
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// The species that have been selected for planting
        /// at each active site.
        /// </summary>
        protected ISiteVar<Dictionary<ISpecies, uint>> PlantingList
        {
            get
            {
                return plantingList;
            }
        }

        //---------------------------------------------------------------------
        /// <summary>
        /// By default, if a form of reproduciton succeeds at a site, it
        /// precludes trying any other forms that haven't been tried yet.
        /// </summary>
        public virtual bool PrecludeRemainingForms
        {
            get {
                return true;
            }
        }

        //---------------------------------------------------------------------

        protected FormOfReproduction()
        {
            int speciesCount = speciesDataset.Count;
            selectedSpecies = Model.Core.Landscape.NewSiteVar<BitArray>();
            plantingList = Model.Core.Landscape.NewSiteVar<Dictionary<ISpecies, uint>>();
            foreach (ActiveSite site in Model.Core.Landscape.ActiveSites) {
                selectedSpecies[site] = new BitArray(speciesCount);
                plantingList[site] = new Dictionary<ISpecies, uint>();
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Are the preconditions for this form of reproduction satisified at
        /// a site for a particular species?
        /// </summary>
        protected abstract bool PreconditionsSatisfied(ISpecies   species, ActiveSite site);

        //---------------------------------------------------------------------

        bool IFormOfReproduction.TryAt(ActiveSite site)
        {
            bool success = false;
            BitArray selectedSpeciesAtSite = selectedSpecies[site];
            Dictionary<ISpecies, uint> plantingListAtSite = plantingList[site];

            for (int index = 0; index < speciesDataset.Count; ++index) {
                if (selectedSpeciesAtSite.Get(index)) {
                    ISpecies species = speciesDataset[index];
                    if (PreconditionsSatisfied(species, site)) {
                        // Temp set propBiomass to 1.0
                        Reproduction.AddNewCohort(species, site,"plant", plantingListAtSite[species]);
                        success = true;
                    }
                }
            }

            return success;
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Clears the list of selected species at a site.
        /// </summary>
        protected void ClearSpeciesAt(ActiveSite site)
        {
            selectedSpecies[site].SetAll(false);
            plantingList[site].Clear();
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Resets the form of reproduction at a site because it will not be
        /// tried.
        /// </summary>
        /// <param name="site">
        /// The site where the form of reproduction will not be tried.
        /// </param>
        void IFormOfReproduction.NotTriedAt(ActiveSite site)
        {
            ClearSpeciesAt(site);
        }
    }
}
