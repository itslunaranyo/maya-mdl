## Maya Export Tool for Quake .MDL

This is the script I use for creating Quake's native model format from within Maya. It evolved slowly from an adaptation of the original MODELGEN command line exporter (which wants Alias .TRI files that no modeling package writes any more), eventually eliminating it entirely and writing straight to .mdl.

This tool requires [Preach's QMDL python library](https://tomeofpreach.wordpress.com/qmdl/).

### Use

Place these files where your Maya install can see them (edit your MAYA.ENV to do so if you haven't already), and import in the Python script window using

```import mdl```

The tool is configured per-model by a text-based script similar to the .qc files used by MODELGEN. Examples are included. "Run" one of these by calling

```mdl.exportScript("pathToFile.txt")```

and a .mdl will be written to the designated place. 

- The +Z axis in Maya becomes +X (east, angle 0) in Quake.
- Animations can be spread across multiple files, referencing a common rig, and read from in sequence using `$file` directives before specifying animation ranges and names with `$anim`. As long as the mesh being animated is an exact match for face count and vertex order, they'll all come together seamlessly. See the Knight example script for more.
- The model can use the old 'squashed base mesh' method of skinning with the `$basemesh` directive, or proper UVs with `$useUVs`.

### Script Format

Export scripts contain commands that begin with a $, one on each line. Any line that doesn't begin with a $command is ignored.

#### Configurations (these go first)

- **$name** - name of model
- **$filedestination** - specify directory to write the .mdl to, final name is $filedestination + $name + .mdl
- **$mesh** - name of mesh object in scene to export. list multiple objects to combine multiple meshes into one.
- **$useUVs** - use the first UV set on the object as the skin coordinates, like a professional
- **$basemesh** - mesh to use for old-school forward-projected UV derivation
- **$skinpadding** - margin to use when using old school $basemesh UVs, in skin pixels
- **$meshnormals** - derive the normals from a different mesh (must have the same vertex counts)
- **$forcenormal** - override the normal of every vert with a particular vector
- **$forwardnode** - locator/object for tracking ai_forward amounts relative to skeleton root
- **$scale** - maya unit and artist mistake compensator
- **$skins** - list of skin names to add to the .mdl (will append $skindir and .bmp)
- **$origin** - specify origin point (monster origins are usually '0 0 24')
- **$flags** - set the flags integer for rocket trails/spinning/etc
- **$projdir** - specify project directory
- **$skindir** - specify directory to look for skins, default is projDir/skins/ otherwise

#### Actions (these go second)

- **$file** - open a maya file. all objects specified by $mesh (and $meshnormals if any) should appear in this file. if $basemesh is used, it only needs to appear in the first $file.
- **$anim** - record vertex positions on all $mesh objects as an animation and add it to the model. format is the anim name, followed by frame numbers or frame ranges (for example, "$anim run 1-10 12" to skip frame 11)
- **$animgroup** - same as $anim but collapse it to a single animgroup in one frame instead of adding all the frames individually. for auto-animating models like torches.
- **$set** - set the value of an attribute in the scene before continuing, format is "$set obj.attr value"

