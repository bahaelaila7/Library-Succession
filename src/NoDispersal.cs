using Landis.Core;
using Landis.SpatialModeling;

namespace Landis.Library.Succession
{
    /// <summary>
    /// Seeding algorithm where no species can seed a neighboring site.
    /// Only the current site is checked.
    /// </summary>
    public static class NoDispersal
    {
        public static void Algorithm(ISpecies species,
                                        ActiveSite site, out bool established, out double seedlingProportion)
        {
            established = Reproduction.SufficientResources(species, site) &&
                   Reproduction.Establish(species, site) &&
                   Reproduction.MaturePresent(species, site);
                   //SiteVars.Cohorts[site].IsMaturePresent(species);
            seedlingProportion = 0;
        }
    }
}
