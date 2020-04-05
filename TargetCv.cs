using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ImputeMissing_NitrogenFixation_Values
{
    public class TargetCv
    {
        public int cutoff;
        public double cv;

        public TargetCv()
        {
            this.cutoff = -1;
            this.cv = double.MaxValue;
        }
    }
}
