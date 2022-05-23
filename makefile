# Build any figure files you find
flist = $(wildcard deconv/figures/figure*.py)
all: $(patsubst deconv/figures/figure%.py, output/figure%.svg, $(flist))

output/figure%.svg: deconv/figures/figure%.py
	@mkdir -p output
	poetry run fbuild $*

test:
	poetry run pytest -s -v -x

clean:
	rm -rf output
