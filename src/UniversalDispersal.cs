using Landis.Core;
using Landis.SpatialModeling;

namespace Landis.Library.Succession
{
    /// <summary>
    /// Seeding algorithm where every species can seed any site; a species does
    /// not even need to be present in any neighboring site.
    /// </summary>
    public static class UniversalDispersal
    {
        public static void Algorithm(ISpecies species,
                                        ActiveSite site, out bool established, out double seedlingProportion)
        {
            established =  Reproduction.SufficientResources(species, site) &&
                   Reproduction.Establish(species, site);
            seedlingProportion = 1;
        }
    }
}
