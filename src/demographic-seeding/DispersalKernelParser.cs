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
using Landis.Core;
using System.Collections.Generic;
//using Landis.Library.Parameters;
using System;

namespace Landis.Library.Succession.DemographicSeeding
{
    /// <summary>
    /// A parser that reads the extension's input and output parameters from
    /// a text file.
    /// </summary>
    public class DispersalKernelParser
        : TextParser<Parameters>
    {
        private ISpeciesDataset speciesDataset;
        private Dictionary<string, int> speciesLineNumbers;

        public static class Names
        {
            public const string LogMessage = "LogMessage";
            public const string ProbabilityTable = "ProbabilityTable";
        }

        //---------------------------------------------------------------------

        static DispersalKernelParser()
        {
            ParsingUtils.RegisterForInputValues();
        }

        //---------------------------------------------------------------------

        public override string LandisDataValue
        {
            get {
                return "Dispersal Kernel";
            }
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public DispersalKernelParser(ISpeciesDataset speciesDataset)
        {
            this.speciesDataset = speciesDataset;
            this.speciesLineNumbers = new Dictionary<string, int>();
        }

        //---------------------------------------------------------------------

        protected override Parameters Parse()
        {
            ReadLandisDataVar();

            Parameters parameters = new Parameters();

            return parameters;
        }

        //---------------------------------------------------------------------

    }
}