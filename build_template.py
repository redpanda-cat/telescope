import os
import argparse
import errno

def merge_files(filename):
    with open(os.path.join("build", "index.html"),"r") as file:
        html = file.read()

    js_txt = open(os.path.join("build", "static/js/main.js"), 'r').read()
    css_txt = open(os.path.join("build", "static/css/main.css"), 'r').read()

    html = html.replace('<script src="/static/js/main.js"></script>', f"<script>{js_txt}</script>")
    html = html.replace('<link href="/static/css/main.css" rel="stylesheet">', f"<style>{css_txt}</style>")

    output_filename = os.path.join("dist", filename, 'template.html')
    if not os.path.exists(os.path.dirname(output_filename)):
        try:
            os.makedirs(os.path.dirname(output_filename))
        except OSError as exc: 
            if exc.errno != errno.EEXIST:
                raise
    with open(output_filename, 'w') as output_file:
        output_file.write(html)


parser = argparse.ArgumentParser(description='Build merged template HTML file')
parser.add_argument('name', help='Name of output file')

args = parser.parse_args()
merge_files(args.name)