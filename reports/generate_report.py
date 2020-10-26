#!/usr/bin/env python
"""
Generate pdf report from the graph

"""

from pathlib import Path
import glob
import re
import markdown
import time
import graph_stat
import argparse
from weasyprint import HTML


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--assembly", help="Graph to report")
    parser.add_argument("-c", "--component", help="Graph to report", nargs="+")
    return parser.parse_args()


def resources_stat(assembly):
    # get the last modified log file
    ak = sorted(Path('logs/construct_graph').glob(pattern=f"asb-{assembly}_*.out"), key=lambda x: x.stat().st_mtime)

    with open(ak[0]) as infile:
        stat_string = (
            "| Statistics | Val(Mb/sec)| Val(GB/mins)|\n"
            "|----|----|----|\n")

        for line in infile:
            if "CPU" in line:
                cpu_time = re.search(r'\d+\.?\d*', line).group()
                stat_string += f"|CPU time| {cpu_time} | {float(cpu_time)/60:.2f}|\n"
            if "Max Memory :" in line:
                max_mem = re.search(r'\d+\.?\d*', line).group()
                stat_string += f"|max mem| {max_mem} | {float(max_mem)/1024:.2f}|\n"
            if "Average Memory :" in line:
                av_mem = re.search(r'\d+\.?\d*', line).group()
                stat_string += f"|average mem| {av_mem} | {float(av_mem)/1024:.2f}|\n"
            if "Run time :" in line:
                run_time = re.search(r'\d+\.?\d*', line).group()
                stat_string += f"|Run time| {run_time} | {float(run_time)/60:.2f}|\n"
    return stat_string


def core_genome_stat(assembly):
    core_string = "### *Core - flexible genome analysis*\n\n\n"
    core_string += "| Pangenome component | Length|Proportion|\n"
    core_string += "| ------------------- | ------|----------|\n"
    core_comp = [x.strip().split() for x in open(f"analysis/core_nonref/{assembly}_core_analysis.tsv").readlines()]
    for param, val in zip(*core_comp):
        core_string += f"| {param} | {int(val):,} | {int(val) *100 / sum( [int(x) for x in core_comp[1][:2]] ):.2f}%\n"
    core_string += "\n\n"
    core_string += f"![Pangenome size as sample increased](analysis/core_nonref/{assembly}_core_flex_sim.png)\n\n\n"
    return core_string


def nonref_stat(assembly):
    nonref_string = "### *Non-reference sequences*\n\n\n"
    nonref_string += "| Non-ref sequences| Node count| Total sequence length|\n"
    nonref_string += "|------------------|-----------|----------------------|\n"
    with open(f"analysis/core_nonref/{assembly}_nonref_analysis.tsv") as infile:
        next(infile)
        for line in infile:
            breed, seq_count, seq_len = line.strip().split()
            nonref_string += f"|{breed}|{int(seq_count):,}|{int(seq_len):,}|\n"
    nonref_string += "\n\n"
    nonref_string += f"![Pangenome nonref sharing count](analysis/core_nonref/{assembly}_nonref_shared_count.png)\n\n\n"
    nonref_string += f"![Pangenome nonref sharing length](analysis/core_nonref/{assembly}_nonref_shared_len.png)\n\n\n"
    return nonref_string


def main():
    args = parse_args()
    assembly = args.assembly
    component = args.component

    intro_string = (f"# **Graph generation report for {assembly}**\n\n\n"
                    f"*Generated on {time.strftime('%a %H:%M %d %B %Y' )}*\n\n"
                    "[TOC]\n"
                    "### *Computational resources*\n")
    resources_string = resources_stat(assembly)
    graph_string = graph_stat.graph_stat_report(assembly, component)
    core_string = core_genome_stat(assembly)
    nonref_string = nonref_stat(assembly)

    all_string = intro_string + resources_string + graph_string + "\n\n" + core_string + nonref_string

    html = markdown.markdown(all_string, extensions=["tables", "toc", "markdown_captions"])
    with open(f"reports/{assembly}_report.html", "w", encoding="utf-8") as outfile:
        outfile.write(html)
        HTML(string=html, base_url=".").write_pdf(f"reports/{assembly}_report.pdf", stylesheets=["report.css"])


if __name__ == "__main__":
    main()
