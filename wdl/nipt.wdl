## Copyright NGeneBio, 2018
## This WDL 
## Requrements/expectations :
## Outputs :
## Cromwell version support
## LICENSING :

# WORKFLOW DEFINITION
workflow NIPTWf {
   File sample_bam
   File? hg19_chrom_size
   File? hg19_bins
   File? hg19_gc
   String? window_size_override

   String default_window_size = select_first([window_size_override, "100000"])
   String bedtools = 'bedtools'
  
   # build bin size
   call CreateBin {
      input:
         hg19_chrom_size = hg19_chrom_size,
         window_size = default_window_size,
         bedtools = bedtools
   }
   
   output{
      File out_bins = CreateBin.out_bins 
   } 
}

#TASK DEFINITIONS
task CreateBin {
   # Command parameters
   File hg19_chrom_size
   String window_size
   String bedtools

   # Runtime paramters
      
   command {
      ${bedtools} makewindows -g ${hg19_chrom_size} -w ${window_size} > hg19_bins 
   }
   output {
      File out_bins = "hg19_bins"
   }
   runtime {
      docker: "biocontainers/bedtools"
      cpu: "1"
      memory: "1GB"
      preemptible: 3
   }
}
