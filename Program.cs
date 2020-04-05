using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using LumenWorks.Framework.IO.Csv;

namespace ImputeMissing_NitrogenFixation_Values
{
    class Program
    {
        static void Main(string[] args)
        {
            //var file = @"H:\Projects\Proteomics\Zymomona\Martein,Julia_NitrogenFixation\R_outputs\EAT_50PercentFilter_NitrogenFixation_DownShiftData_ZM4tag.csv";
            var file = @"H:\Projects\Proteomics\Zymomona\Martein,Julia_NitrogenFixation\R_outputs\EAT_50PercentFilter_NitrogenFixation_UpShiftData_ZM4tag.csv";
            var reader = new CsvReader(new StreamReader(file), true, ',');

            var headers = reader.GetFieldHeaders();
            var replicates = new List<Replicate>();

            // read lfq spreadsheet. Parse data into replicate objects.
            // replicate data objects hold all lfq values for a sample and metadata
            // concerning replicate/condition names.
            foreach (var header in headers)
            {
                Console.WriteLine(header + " : " + reader.GetFieldIndex(header));

                // h3k spreadsheet logic
                //if (header.Contains("LFQ intensity"))
                if (header.Contains("0") || header.Contains("5") || header.Contains("15") || header.Contains("30") || header.Contains("60") || header.Contains("120") || header.Contains("NH4") || header.Contains("N2")) //|| header.Contains("N2-120"))
                {
                    replicates.Add(new Replicate(header, reader.GetFieldIndex(header)));
                    //Console.WriteLine(replicates);
                }
            }

            while (reader.ReadNextRecord())
            {
                foreach (var replicate in replicates)
                {
                    double lfqValue = -1;

                    try
                    {
                        if (double.TryParse(reader[replicate.readerIndex], out lfqValue))
                        {
                            replicate.lfqValues.Add(new LfqValues(Math.Log(lfqValue, 2), false, reader[0]));
                        }
                        else
                        {
                            replicate.lfqValues.Add(new LfqValues(0, true, reader[0]));
                        }
                    }
                    catch (MissingFieldCsvException e)
                    {
                        replicate.lfqValues.Add(new LfqValues(0, true, reader[0]));
                    }
                }
            }

            var uniqueConditions = replicates.Select(replicate => replicate.conditionName).ToList().Distinct();
            var conditions = new List<Condition>();

            foreach (var condition in uniqueConditions)
            {
                var newCondition = new Condition(condition);
                newCondition.replicates = replicates.Where(replicate => replicate.conditionName.Equals(condition)).ToList();

                conditions.Add(newCondition);
                //Console.WriteLine(condition);
            }

            //conditions.Last().RunImputationAlgorithm();


            foreach (var condition in conditions)
            {
                condition.RunImputationAlgorithm();
            }

            // H3K Spreadsheet Logic
            // now write out results file
            // two files, one summarizing the cutoffs for each condition
            /*
            var t = "";

            var outputPath = @"P:\IJM_H3K\Imputed Data\ImputationCutoffs_11%.csv";
            var writer = new StreamWriter(outputPath);

            // write headers
            writer.WriteLine("Condition, Cutoff");
            foreach (var condition in conditions)
            {
                writer.WriteLine("{0},{1}", condition.conditionName, condition.imputationCutoff);
            }

            writer.Close();
            writer.Dispose();
            */

            // now write out results file
            // two files, one summarizing the cutoffs for each condition
            var outputPath = @"H:\Projects\Proteomics\Zymomona\Martein,Julia_NitrogenFixation\R_outputs\NitrogenFixation_Upshift_ZM4tag_imputation_cutoffs.csv";
            var writer = new StreamWriter(outputPath);

            // write headers
            writer.WriteLine("Condition, Cutoff %");
            foreach (var condition in conditions)
            {

                writer.WriteLine("{0},{1}", condition.conditionName, condition.imputationCutoff / 10);
            }

            writer.Close();
            writer.Dispose();

            // write out file with imputed values
            outputPath = @"H:\Projects\Proteomics\Zymomona\Martein,Julia_NitrogenFixation\R_outputs\NitrogenFixation_UpshiftImputed_ZM4tag_proteinGroups.csv";
            writer = new StreamWriter(outputPath);

            // write headers
            WriteHeaders(writer, conditions);

            // write data
            WriteData(writer, conditions);

            writer.Close();
            writer.Dispose();
        }

        public static void WriteHeaders(StreamWriter writer, List<Condition> conditions)
        {
            writer.Write("Majority protein ID");

            foreach (var condition in conditions)
            {
                foreach (var replicate in condition.replicates)
                {
                    writer.Write(",LFQ intensity {0}_{1}", condition.conditionName, replicate.replicateName);
                }
            }

            writer.WriteLine();
        }

        public static void WriteData(StreamWriter writer, List<Condition> conditions)
        {
            for (var i = 0; i < Condition.biomolecules.Count; i++)
            {
                WriteData(writer, conditions, i);
            }
        }

        private static void WriteData(StreamWriter writer, List<Condition> conditions, int index)
        {
            writer.Write(Condition.biomolecules[index]);

            foreach (var condition in conditions)
            {
                foreach (var replicate in condition.replicates)
                {
                    writer.Write(",{0}", Math.Pow(2, replicate.lfqValues[index].lfq));
                }
            }

            writer.WriteLine();
        }
    }
}
