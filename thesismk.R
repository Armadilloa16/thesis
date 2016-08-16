
library(knitr)

main_file = 'whole_thesis'
perl_dir = 'C:\\Perl\\bin\\'
tex_dir = 'C:\\Program Files (x86)\\MiKTeX 2.9\\'
# tex_dir = 'C:\\Program Files\\MiKTeX 2.9\\'

setwd("texput")
unlink(paste(main_file,"*",sep=""))
unlink('knitr.sty')
unlink('./logs/*.log')
knit(paste('../',main_file,'.Rnw',sep=""),
     output=paste(main_file,'.tex',sep="")
)
# system2('pdflatex',
#         args=main_file,
#         stdout = paste('./logs/',main_file,'_pdflatex_1_stdout.log',sep=""),
#         stderr = paste('./logs/',main_file,'_pdflatex_1_stderr.log',sep=""),
#         wait = TRUE
# )
# system2('bibtex',
#         args=main_file,
#         stdout = paste('./logs/',main_file,'_bibtex_stdout.log',sep=""),
#         stderr = paste('./logs/',main_file,'_bibtex_stderr.log',sep=""),
#         wait = TRUE
# )
# system2('pdflatex',
#         args=main_file,
#         stdout = paste('./logs/',main_file,'_pdflatex_2_stdout.log',sep=""),
#         stderr = paste('./logs/',main_file,'_pdflatex_2_stderr.log',sep=""),
#         wait = TRUE
# )
# system2(paste(perl_dir,'perl',sep=""),
#         args=paste('"',tex_dir,'scripts\\glossaries\\makeglossaries','" ',main_file,sep=""),
#         stdout = paste('./logs/',main_file,'_makeglossaries_stdout.log',sep=""),
#         stderr = paste('./logs/',main_file,'_makeglossaries_stderr.log',sep=""),
#         wait = TRUE
# )
# system2('pdflatex',
#         args=main_file,
#         stdout = paste('./logs/',main_file,'_pdflatex_3_stdout.log',sep=""),
#         stderr = paste('./logs/',main_file,'_pdflatex_3_stderr.log',sep=""),
#         wait = TRUE
# )
setwd("..")

