import MDAnalysis as mda
import click

def set_chain_id(pdb_path: str, output_path: str,chain_id: str):
    lines = []
    with open(pdb_path, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                lines.append(line[:21] + chain_id + line[22:])
            else:
                lines.append(line)

    with open(output_path, 'w') as f:
        f.writelines(lines)

@click.command()
@click.option('--pdb_path', '-i', help='Input pdb file path')
@click.option('--output_path', '-o', help='Output pdb file path')
@click.option('--chain_id', '-c', help='Chain ID')
def set_chain_id_to_command(pdb_path: str, output_path: str,chain_id: str):
    set_chain_id(pdb_path, output_path,chain_id)
