TARGET=kernel_generator_spec
SRC=$(TARGET).tex
PDF=$(TARGET).pdf

all:	$(PDF)

$(PDF):	$(SRC)
	pdflatex $(SRC)
clean:
	rm $(TARGET).aux $(TARGET).log $(TARGET).toc
