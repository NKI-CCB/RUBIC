# Copyright 2015 Netherlands Cancer Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


focal.events.to.tsv <- function(focal.events, file.name='') {
  if (isempty(file.name)) {
    file.name <- stdout()
  }
  
  fd <- file(file.name, "wb")
  lines <- vapply(focal.events, function(event) {
    paste(event$chromosome, event$loc.start, event$loc.end,
          formatC(event$l$q, format="e", digits=5),
          formatC(event$r$q, format="e", digits=5),
          paste0(event$gene.symbols, collapse=","), sep='\t')
  }, character(1))
  lines <- c(paste('Chromosome', 'Start', 'End', 'Left break (-log10(qValue))',
                   'Right break (-log10(qValue))', 'Gene symb', sep='\t'),
             lines)
  writeLines(lines, con=fd, sep='\n')
  flush(fd)
  close(fd)
}