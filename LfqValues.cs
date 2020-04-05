using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ImputeMissing_NitrogenFixation_Values
{
    public class LfqValues
    {
        public string bioMolecule;
        public double lfq;
        public bool isImputed;

        public LfqValues(double lfq, bool isImputed, string bioMolecule)
        {
            this.lfq = lfq;
            this.isImputed = isImputed;
            this.bioMolecule = bioMolecule;
        }

        public override string ToString()
        {
            return this.lfq.ToString();
        }
    }
}
