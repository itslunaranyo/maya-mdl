## Maya Export Tool for Quake .MDL

This is the script I use for creating Quake's native model format from within Maya. It evolved slowly from an adaptation of the original MODELGEN command line exporter (which wants Alias .TRI files that no modeling package writes any more), eventually eliminating it entirely and writing straight to .mdl.

This tool requires [Preach's QMDL python library](https://tomeofpreach.wordpress.com/qmdl/).

### Use

Place these files where your Maya install can see them (edit your MAYA.ENV to do so if you haven't already), and import in the Python script window using

```import mdl```

The tool is configured per-model by a text-based script similar to the .qc files used by MODELGEN. Examples are included. "Run" one of these by calling

```mdl.exportScript("pathToFile.txt")```

and a .mdl will be written to the designated place. 

### Notes

- The +Z axis in Maya becomes +X (east, angle 0) in Quake.
- Animations can be spread across multiple files, referencing a common rig, and read from in sequence using `$file` directives before specifying animation ranges and names with `$anim`. As long as the mesh being animated is an exact match for face count and vertex order, they'll all come together seamlessly. See the Knight example script for more.
- The model can use the old 'squashed base mesh' method of skinning with the `$basemesh` directive, or proper UVs with `$useUVs`.