#
# $Id: symmetry.tcl,v 1.15 2013/04/15 17:45:51 johns Exp $
#
#VMD is Copyright © 1995-2016 Theoretical and Computational Biophysics Group and the
#Board of Trustees of the University of Illinois
#Portions of this code are copyright © 1997-1998 Andrew Dalke.
package provide symmetrytool 1.3

#puts "*** symmetrytool ***"
#puts "*** developed by VMD ***"
#puts "*** modified by Jinglin Mu ***"
#puts "*** watfros@sdut.edu.cn ***"


namespace eval ::Symmetry:: {

   proc initialize {} {
      variable tol_inversion
      variable tol_rotation
      variable tol_reflection
      variable canvaso -1
      variable canvast -1
      
      variable selectcolor  lightsteelblue
      variable fixedfont    {Courier 9};
      variable bigfont      {Helvetica 12}
      variable bigboldfont  {Helvetica 12 bold}
      variable axeslistformat [format " #  %5s" order]

      variable showunique  0
      variable showideal   0
      variable showorient  0
      variable showinertia 0
      variable showcenter  1
      variable showrotaxes 1
      variable showrraxes  1
      variable showplanes  1
      variable checkbonds  1

      variable lastselectaxis 0 
      variable lastselectplaneo 0 
      variable lastselectplanet 0 
      variable lastselectrraxis 0 
      variable axislist {}
      variable axisgid {}
      variable rraxislist {}
      variable rraxisgid {}
      variable planelist {}
      variable planegid {}
      variable centerinfo {}
      variable centerid {}
      variable inertialist {}
      variable inertiagido {}
      variable inertiagidi {}
      variable inversion 0 

      variable repname    {}
      variable uniquerepo {}
      variable uniquerepi {}
      variable origmol    -1
      variable idealmol   -1
      variable origmolname none

      variable tol      0.05
      variable elements {}
      variable orient   {}
      variable unique   {}
      variable ideal    {}
      variable rmsd     ---
      variable com      {}
      variable inertia  {}
      variable seltext  all
      variable sel      noselection

      variable pointgroup -none-
   }

   initialize
}

##
## Main routine
## Create the window and initialize data structures
##
proc symmetrytool_tk {} {
  Symmetry::create_gui
  return .symmetry
}

proc ::Symmetry::create_gui {args} {
   initialize

   set sel [::util::getargs $args "-sel" {}]
   if {[llength $sel]} {
      variable seltext [$sel text]
      variable origmol [$sel molid]
      make_selection
   } else {
      #variable origmol [molinfo top]
      # Set to first molecule that has coordinates
      foreach mol [molinfo list] {
         if {[molinfo $mol get numframes]>0 &&
             $mol!=$::Symmetry::idealmol} {
            variable origmolname "$mol [molinfo $mol get name]"
            variable origmol $mol
            break
         }
      }
   }

   if {[winfo exists .symmetry]} {
      return
   }

   set w [toplevel ".symmetry"]
   wm title $w "SymmetryTool"
   wm resizable $w 1 1

   wm protocol .symmetry WM_DELETE_WINDOW {
      # Delete the selection
      if {![catch {$::Symmetry::sel num}]} {
         $::Symmetry::sel delete
      }
      # Delete the drawing canvases
      if {$::Symmetry::canvaso>=0 && 
          [lsearch [molinfo list] $::Symmetry::canvaso]>=0} {
         mol delete $::Symmetry::canvaso
         set ::Symmetry::canvaso -1
      }
      if {$::Symmetry::canvast>=0 && 
          [lsearch [molinfo list] $::Symmetry::canvast]>=0} {
         mol delete $::Symmetry::canvast
         set ::Symmetry::canvast -1
      }
      if {$::Symmetry::idealmol>=0 && 
          [lsearch [molinfo list] $::Symmetry::idealmol]>=0} {
         mol delete $::Symmetry::idealmol
         set ::Symmetry::idealmol -1
      }
      set molid $::Symmetry::origmol
      if {$molid>=0} {
         set repindex [mol repindex $molid $::Symmetry::repname]
         if {$repindex>=0} {
            mol delrep $repindex $molid
         }
         set repindex [mol repindex $molid $::Symmetry::uniquerepo]
         if {$repindex>=0} {
            mol delrep $repindex $molid
         }
         foreach gid $::Symmetry::inertiagido {
            graphics $::Symmetry::origmol delete $gid
         }
         set numreps [molinfo $molid get numreps]
         for {set i 0} {$i<$numreps} {incr i} {
            mol showrep $molid $i 1
         }
         if {![catch [list $::Symmetry::sel num]]} {
            $::Symmetry::sel set {x y z} $::Symmetry::origcoor
         }
         mol on $molid
      }
      destroy .symmetry
   }

   # System molecule
   frame $w.head
   set f $w.head
   label $f.mollabel -text "Molecule:" -anchor w
   menubutton $f.inputmol -relief raised -bd 2 -direction flush \
       -textvariable ::Symmetry::origmolname -width 30 \
       -menu $f.inputmol.menu
   menu $f.inputmol.menu -tearoff no
   grid $f.mollabel  -in $f -row 1 -column 0 -columnspan 1 -sticky w -pady 2 -padx 2
   grid $f.inputmol  -in $f -row 1 -column 1 -columnspan 1 -sticky w -pady 2  
   pack $f -anchor w 

   frame $w.sel
   label $w.sel.label -text "Atom selection: "
   entry $w.sel.entry -textvariable ::Symmetry::seltext
   pack $w.sel.label -side left
   pack $w.sel.entry -side left -expand 1 -fill both

   variable bigfont
   variable bigboldfont
   frame $w.guess
   label $w.guess.label  -font $bigfont -text "Point group:"
   label $w.guess.entry  -font $bigboldfont -width 6 -textvariable ::Symmetry::pointgroup
   button $w.guess.button -text "Guess-Symmetry" -command ::Symmetry::make_guess
   label $w.guess.tollabel -text "Tolerance: "
   entry $w.guess.tolentry -width 5 -textvariable ::Symmetry::tol

   pack $w.guess.label -side left
   pack $w.guess.entry $w.guess.button -side left -padx 10
   pack $w.guess.tollabel $w.guess.tolentry -side left

   checkbutton $w.bonds -text "Consider bond order and orientation in guess" \
      -variable ::Symmetry::checkbonds

   frame $w.summary
   label $w.summary.label -text "Symmetry elements: "
   label $w.summary.entry -textvariable ::Symmetry::elements
   pack $w.summary.label $w.summary.entry -anchor w -side left

   frame $w.rmsd
   label $w.rmsd.label -text "RMSD between original and idealized structure:"
   label $w.rmsd.entry -textvariable ::Symmetry::rmsd
   pack  $w.rmsd.label $w.rmsd.entry  -side left

   frame $w.m
   frame $w.m.check
   #checkbutton $w.m.check.unique -text "Show unique atoms" -variable ::Symmetry::showunique \
      -command ::Symmetry::toggle_unique
   #pack $w.m.check.unique -anchor w

   #checkbutton $w.m.check.ideal -text "Show idealized coordinates" \
      -variable ::Symmetry::showideal -command ::Symmetry::toggle_ideal
   #pack $w.m.check.ideal -anchor w

   #checkbutton $w.m.check.orient -text "Standard orientation" \
      -variable ::Symmetry::showorient -command ::Symmetry::toggle_orient
   #pack $w.m.check.orient -anchor w

   #checkbutton $w.m.check.princ -text "Show principle axes of inertia" \
      -variable ::Symmetry::showinertia -command ::Symmetry::toggle_inertia
   #pack $w.m.check.princ -anchor w


   variable selectcolor
   variable fixedfont

   #labelframe $w.m.inert -bd 2 -relief ridge -text "Principle axes of inertia" -padx 2m -pady 2m
   #label $w.m.inert.format -font $fixedfont -text "\#  eigenval  unique" \
   #   -relief flat -bd 2 -justify left;
   #listbox $w.m.inert.list -font $fixedfont \
      -width 20 -height 3 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::Symmetry::inertialist

   #pack $w.m.inert.format  -anchor w
   #pack $w.m.inert.list -side left -anchor w 

   #pack $w.m.check -side left -anchor w 
   #pack $w.m.inert -side right -anchor w -padx 1m
  
   frame $w.ele

   # Rotary axes
   labelframe $w.ele.axes -bd 2 -relief ridge -text "Rotary axes" -padx 2m -pady 2m

   checkbutton $w.ele.axes.show -text "Show" \
      -variable ::Symmetry::showrotaxes -command ::Symmetry::draw_symmetry_elements
   pack $w.ele.axes.show -anchor w

   label $w.ele.axes.format -font $fixedfont -text " \# order type" \
      -relief flat -bd 2 -justify left;

   frame $w.ele.axes.list
   scrollbar $w.ele.axes.list.scroll -command "$w.ele.axes.list.list yview"
   listbox $w.ele.axes.list.list -yscroll "$w.ele.axes.list.scroll set" -font $fixedfont \
      -width 20 -height 12 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::Symmetry::axislist

   pack $w.ele.axes.list.list    -side left -fill both -expand 1
   pack $w.ele.axes.list.scroll  -side left -fill y    -expand 0

   pack $w.ele.axes.format  -anchor w
   pack $w.ele.axes.list    -expand 1 -fill both

   # Rotary reflections
   labelframe $w.ele.rraxes -bd 2 -relief ridge -text "Rotary reflections" -padx 2m -pady 2m

   checkbutton $w.ele.rraxes.show -text "Show" \
      -variable ::Symmetry::showrraxes -command ::Symmetry::draw_symmetry_elements
   pack $w.ele.rraxes.show -anchor w

   label $w.ele.rraxes.format -font $fixedfont -text " \# order type" \
      -relief flat -bd 2 -justify left;

   frame $w.ele.rraxes.list
   scrollbar $w.ele.rraxes.list.scroll -command "$w.ele.rraxes.list.list yview"
   listbox $w.ele.rraxes.list.list -yscroll "$w.ele.rraxes.list.scroll set" -font $fixedfont \
      -width 17 -height 12 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::Symmetry::rraxislist

   pack $w.ele.rraxes.list.list    -side left -fill both -expand 1
   pack $w.ele.rraxes.list.scroll  -side left -fill y    -expand 0

   pack $w.ele.rraxes.format  -anchor w
   pack $w.ele.rraxes.list    -expand 1 -fill both

   # Planes
   labelframe $w.ele.planes -bd 2 -relief ridge -text "Mirror planes" -padx 2m -pady 2m

   checkbutton $w.ele.planes.show -text "Show" \
      -variable ::Symmetry::showplanes -command ::Symmetry::draw_symmetry_elements
   pack $w.ele.planes.show -anchor w

   label $w.ele.planes.format -font $fixedfont -text " \# type" \
      -relief flat -bd 2 -justify left;

   frame $w.ele.planes.list
   scrollbar $w.ele.planes.list.scroll -command "$w.ele.planes.list.list yview"
   listbox $w.ele.planes.list.list -yscroll "$w.ele.planes.list.scroll set" -font $fixedfont \
      -width 18 -height 12 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::Symmetry::planelist

   pack $w.ele.planes.list.list    -side left -fill both -expand 1
   pack $w.ele.planes.list.scroll  -side left -fill y    -expand 0

   pack $w.ele.planes.format  -anchor w
   pack $w.ele.planes.list    -expand 1 -fill both

   labelframe $w.ele.center -bd 2 -relief ridge -text "Centre of symmetry" -padx 2m -pady 2m

   #pack $w.ele.axes $w.ele.rraxes $w.ele.planes -padx 1m -pady 1m -side left -expand 1 -fill both
   # 对称中心 labelframe 内的组件
   checkbutton $w.ele.center.show -text "Show" \
      -variable ::Symmetry::showcenter -command ::Symmetry::draw_symmetry_elements
   pack $w.ele.center.show -anchor w

   label $w.ele.center.format -font $fixedfont -text " \# inversion" \
      -relief flat -bd 2 -justify left;

   frame $w.ele.center.list
   scrollbar $w.ele.center.list.scroll -command "$w.ele.center.list.list yview"
   listbox $w.ele.center.list.list -yscroll "$w.ele.center.list.scroll set" -font $fixedfont \
      -width 18 -height 12 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::Symmetry::centerinfo

   pack $w.ele.center.list.list    -side left -fill both -expand 1
   pack $w.ele.center.list.scroll  -side left -fill y    -expand 0

   pack $w.ele.center.format  -anchor w
   pack $w.ele.center.list    -expand 1 -fill both

   # 使用 grid 进行 2x2 布局
   grid $w.ele.axes -row 0 -column 0 -padx 1m -pady 1m -sticky "nsew"
   grid $w.ele.rraxes -row 0 -column 1 -padx 1m -pady 1m -sticky "nsew"
   grid $w.ele.planes -row 1 -column 0 -padx 1m -pady 1m -sticky "nsew"
   grid $w.ele.center -row 1 -column 1 -padx 1m -pady 1m -sticky "nsew"

   # 设置网格的权重，使组件能均匀扩展
   grid rowconfigure $w.ele 0 -weight 1
   grid rowconfigure $w.ele 1 -weight 1
   grid columnconfigure $w.ele 0 -weight 1
   grid columnconfigure $w.ele 1 -weight 1

   pack $w.sel -padx 1m -pady 1m -expand 1 -fill both
   pack $w.guess $w.rmsd $w.summary -padx 1m -pady 1m -anchor w
   pack $w.m $w.ele -padx 1m -pady 1m -expand 1 -fill both
 
   # This will be executed when items are selected:   
   bind $w.ele.axes.list.list <<ListboxSelect>> {
      ::Symmetry::select_axis [.symmetry.ele.axes.list.list curselection]
   }

   bind $w.ele.rraxes.list.list <<ListboxSelect>> {
      ::Symmetry::select_rraxis [.symmetry.ele.rraxes.list.list curselection]
   }

   bind $w.ele.planes.list.list <<ListboxSelect>> {
      ::Symmetry::select_plane [.symmetry.ele.planes.list.list curselection]
   }

   bind $w.ele.center.list.list <<ListboxSelect>> {
      ::Symmetry::select_center [.symmetry.ele.center.list.list curselection]
   }
   # Update panes
   update_mol_list

   trace add variable ::vmd_trajectory_read   write ::Symmetry::update_mol_list
   trace add variable ::vmd_molecule          write ::Symmetry::update_mol_list
   trace add variable ::Symmetry::origmolname write ::Symmetry::update_mol_menu
   trace add variable ::Symmetry::origmol     write ::Symmetry::update_mol
}


proc Symmetry::update_mol_menu {args} {
   set m [lindex $::Symmetry::origmolname 0]
   if {$m=="none"} {
      set ::Symmetry::origmol -1
   } else {
      set ::Symmetry::origmol $m
   }
}

# traced command to autoupdate menus when number of mols is changed
proc Symmetry::update_mol_list {args} {
   set mollist [molinfo list]
   set some_mols_have_coords 0

   variable idealmol
   if {[lsearch $mollist $::Symmetry::idealmol]<0} {
      set idealmol -1
   }
   if {[lsearch $mollist $::Symmetry::canvaso]<0} {
      set ::Symmetry::canvaso -1
   }
   if {[lsearch $mollist $::Symmetry::canvast]<0} {
      set ::Symmetry::canvast -1
   }

   set f .symmetry.head
   set menu $f.inputmol.menu
   $menu delete 0 last
   $f.inputmol configure -state disabled

   variable origmol
   variable origmolname
   if {($origmol<0 && [llength $mollist]) ||
       [lsearch $mollist $origmol]<0} {
      # Set to first molecule that has coordinates
      foreach mol $mollist {
         if {[molinfo $mol get numframes]>0 &&
             $mol!=$idealmol} {
            set origmolname "$mol [molinfo $mol get name]"
            break
         }
      }
   }

   foreach mol $mollist {
      $menu add radiobutton -label "$mol [molinfo $mol get name]" \
         -variable ::Symmetry::origmolname
      if {[molinfo $mol get numframes] > 0 &&
          $mol!=$idealmol} {
         set some_mols_have_coords 1
      } else {
         $menu entryconfigure [expr $mol + 0] -state disabled
      }
   }

   if {$some_mols_have_coords} {
      $f.inputmol configure -state normal
   } elseif {![llength $mollist]} {
      set origmolname none
   }
}


proc ::Symmetry::update_mol {args} {
   variable canvaso
   variable canvast
   variable lastselectplaneo {}
   variable lastselectplanet {}
   variable lastselectaxis   {}
   variable lastselectrraxis {}
   if {$canvaso>=0} { catch [list graphics $canvaso delete all] }
   if {$canvast>=0} { catch [list graphics $canvast delete all] }
   variable axisgid   {}
   variable planegid  {}
   variable rraxisgid {}
   variable axes {}
   variable axislist {}
   variable rraxes {}
   variable rraxislist {}
   variable planes {}
   variable planelist {}
   variable centerinfo {}
   variable centerid {}
   variable inertialist {}
   variable inertiagido {}
   variable inertiagidi {}
   
   variable repname    {}
   variable uniquerepo {}
   variable uniquerepi {}

   variable pointgroup -none-
   variable elements {}
   variable orient   {}
   variable unique   {}
   variable ideal    {}
   variable rmsd     ---
   variable com      {}
   variable inertia  {}

   catch [list $::Symmetry::sel delete]
   set ::Symmetry::sel noselection
}


proc ::Symmetry::make_guess {} {
   if {![llength $::Symmetry::seltext]} { return }

   ::Symmetry::make_selection
   if {![catch [list $::Symmetry::sel num]] && [$::Symmetry::sel num]} {
      vmdcon "Guessing symmetry"
      ::Symmetry::guess $::Symmetry::sel $::Symmetry::tol
   }
}

proc ::Symmetry::make_selection {} {
   if {$::Symmetry::origmol<0} { return }
   if {![llength $::Symmetry::seltext]} { return }

   uplevel \#0 {
      #puts $::Symmetry::seltext
      set oldsel $::Symmetry::sel

      set ::Symmetry::sel [atomselect $::Symmetry::origmol $::Symmetry::seltext]
      $::Symmetry::sel global
      #puts $::Symmetry::sel
      
      # Delete the selection
      if {![catch {$oldsel num}]} {
         $oldsel delete
      }
   }

   if {[catch {$::Symmetry::sel num}]} { return }

   variable repname
   variable seltext
   variable origmol
   if {[mol repindex $origmol $repname]<0} {
      mol selection $seltext
      mol representation Licorice 0.1 10 10
      mol addrep $origmol
      variable repname [mol repname $origmol [expr {[molinfo $origmol get numreps]-1}]]
   } else {
      mol modselect [expr {[mol repindex $origmol $repname]}] $origmol $seltext
   }

   set numreps [molinfo $origmol get numreps]
   for {set i 0} {$i<$numreps} {incr i} {
      mol showrep $origmol $i 0
   }
   mol showrep $origmol [mol repindex $origmol $repname] 1

   variable idealmol
   if {$idealmol>0} {
      mol modselect 0 $idealmol $seltext      
   }
}


# 选中对称中心的处理过程
proc ::Symmetry::select_center {item} {
   variable centerid
   variable canvaso
   if {![llength $item] || $canvaso<0} { return }
   #puts "Selected item $item"

   variable lastselectcenter
   # 先将对称中心恢复默认颜色
   graphics $canvaso replace $lastselectcenter
   graphics $canvaso color purple ;# 假设默认颜色是紫色

   set lastselectcenter [lindex $centerid $item 0]
   # 将选中的对称中心颜色变为黄色以高亮显示
   graphics $canvaso replace $lastselectcenter
   graphics $canvaso color yellow

   # 更新显示的详细信息，这里简单示例打印信息
   puts "Selected center"
   set allatom [atomselect top all]
   set natom [$allatom num]
   #每隔0.1秒执行一个原子平移动画，执行10次
   for {set i 1} {$i <= $natom} {incr i} {
   array set atom_ ""
   array set vscale_ ""
   array set vatom_ ""
   set atom_($i) [atomselect top "serial $i"]
   set vatom_($i) [lindex [$atom_($i) get {x y z}] 0]
   #puts "vatom_($i)"
   }
   
   for {set j 1} {$j <= 10} {incr j} {
   for {set i 1} {$i <= $natom} {incr i} {
   $atom_($i) moveby [vecscale 0.2 [vecinvert $vatom_($i)]]
   }
   display update
   after 100
   }

}

proc ::Symmetry::select_axis {item} {
   variable canvaso
   variable axes
   variable sel
   if {![llength $item] || $canvaso<0} { return }
   
   variable axisgid
   variable lastselectaxis
   graphics $canvaso replace $lastselectaxis
   graphics $canvaso color orange3

   set lastselectaxis [lindex $axisgid $item 0]

   graphics $canvaso replace $lastselectaxis
   graphics $canvaso color orange2
   
   
   #每隔0.1秒执行一个旋转动画，执行10次
   set axis_i [lindex $axes $item]
   #puts "Selected axis [lindex $axis_i 0]"
   set output [lindex $axis_i 1]
   set output [expr {36.00 / $output }]
   set rdeg $output
   set allatom [atomselect top all]
   set natom [$sel num]
   for {set i 1} {$i <= 10} {incr i} {
   display update
   $allatom move [trans axis [lindex $axis_i 0] $rdeg deg]
   after 100
   }
}

proc ::Symmetry::select_rraxis {item} {
   variable canvaso
   if {![llength $item] || $canvaso<0} { return }

   puts "Selected item $item"
   variable rraxisgid
   variable lastselectrraxis
   variable rraxes
   variable sel
   
   graphics $canvaso replace $lastselectrraxis
   graphics $canvaso color lime

   set lastselectrraxis [lindex $rraxisgid $item 0]

   graphics $canvaso replace $lastselectrraxis
   graphics $canvaso color cyan2
   
   #每隔0.1秒执行一个旋转动画，执行10次
   #puts "Selected item [lindex $axes $item]"
   set rraxis_i [lindex $rraxes $item]
   #puts "Selected axis [lindex $axis_i 0]"
   set output [lindex $rraxis_i 1]
   set output [expr {36.00 / $output }]
   set rdeg $output
   #这一句实际上描述了初始载入分子的数目
   #set natom [$sel num]
   #如何想通过使用某个分子的对称元素，对另外一个分子变形，那么必须使用下面的语句进行处理
   set allatom [atomselect top all]
   set natom [$allatom num]
   #puts $natom
   
   for {set i 1} {$i <= 10} {incr i} {
   display update
   $allatom move [trans axis [lindex $rraxis_i 0] $rdeg deg]
   after 100
   }
   
   #每隔0.1秒执行一个原子平移动画，执行10次
   for {set i 1} {$i <= $natom} {incr i} {
   array set atom_ ""
   array set vscale_ ""
   array set vatom_ ""
   set atom_($i) [atomselect top "serial $i"]
   set vatom_($i) [lindex [$atom_($i) get {x y z}] 0]
   set vscale_($i) [vecdot [lindex $rraxis_i 0] $vatom_($i)]
   set vscale_($i) [expr {0.2*$vscale_($i)}]
   }
   
   for {set j 1} {$j <= 10} {incr j} {
   for {set i 1} {$i <= $natom} {incr i} {
   $atom_($i) moveby [vecscale $vscale_($i) [vecinvert [lindex $rraxis_i 0]]]
   }
   display update
   after 100
   }
}

proc ::Symmetry::select_plane {item} {
   variable canvaso
   variable canvast
   if {![llength $item] || $canvaso<0 || $canvast<0} { return }

   puts "Selected item $item"
   variable planegid
   variable lastselectplaneo
   variable lastselectplanet
   variable planes
   variable sel
   
   graphics $canvaso replace $lastselectplaneo
   graphics $canvaso color pink
   graphics $canvast replace $lastselectplanet
   graphics $canvast color pink

   set lastselectplaneo [lindex $planegid $item 0 0]
   set lastselectplanet [lindex $planegid $item 1 0]

   graphics $canvaso replace $lastselectplaneo
   graphics $canvaso color magenta
   graphics $canvast replace $lastselectplanet
   graphics $canvast color magenta
   
   #每隔0.1秒执行一个原子平移动画，执行10次
   #这一句实际上描述了初始载入分子的原子数目
   #set natom [$sel num]
   #如何想通过使用某个分子的对称元素，对另外一个分子变形，那么必须使用下面的语句进行处理
   set allatom [atomselect top all]
   set natom [$allatom num]
   #puts $natom
   
   set planes_i [lindex $planes $item]
   for {set i 1} {$i <= $natom} {incr i} {
   array set atom_ ""
   array set vscale_ ""
   array set vatom_ ""
   set atom_($i) [atomselect top "serial $i"]
   set vatom_($i) [lindex [$atom_($i) get {x y z}] 0]
   set vscale_($i) [vecdot [lindex $planes_i 0] $vatom_($i)]
   set vscale_($i) [expr {0.2*$vscale_($i)}]
   }
   
   for {set j 1} {$j <= 10} {incr j} {
   for {set i 1} {$i <= $natom} {incr i} {
   $atom_($i) moveby [vecscale $vscale_($i) [vecinvert [lindex $planes_i 0]]]
   }
   display update
   after 100
   }
   
}

proc ::Symmetry::toggle_ideal {} {
   variable showideal
   variable idealmol
   variable origmol
   variable sel
   if {[catch {$sel num}] || $idealmol<0 || $origmol<0} { return }

   if {$showideal} {
      mol on  $idealmol 
      mol off $origmol  
   } else {
      mol off $idealmol 
      mol on  $origmol  
   }
   toggle_unique
}

proc ::Symmetry::toggle_unique {} {
   variable showunique
   variable showideal
   variable idealmol
   variable origmol
   variable uniquerepo
   variable uniquerepi
   if {$origmol<0} { return }

   set mol $origmol
   set uniquerep $uniquerepo
   if {$showideal} {
      set mol $idealmol
      set uniquerep $uniquerepi
   }

   if {![llength $uniquerep]} {
      set irep [molinfo $mol get numreps]
      variable unique
      mol selection "index $unique"
      mol representation Licorice 0.3 10 10 
      mol addrep $mol
      set uniquerep [mol repname $mol $irep]
      if {$showideal} {
         set uniquerepi $uniquerep
      } else {
         set uniquerepo $uniquerep
      }
   }

   if {$showunique} {
      mol showrep $mol [mol repindex $mol $uniquerep] 1
   } elseif {[llength $uniquerep]} {
      mol showrep $mol [mol repindex $mol $uniquerep] 0
   }
}

proc ::Symmetry::toggle_orient {} {
   variable orient
   if {![llength $orient]} { return }
   variable origmol
   variable sel
   if {$origmol<0 || [catch {$sel num}]} { return }

   variable showorient
   variable seltext
   variable idealmol

   set all [atomselect $idealmol all]
   if {$showorient} {
      $sel move $orient
      $all move $orient
   } else {
      variable idealcoor
      $all set {x y z} $idealcoor
      variable origcoor
      $sel set {x y z} $origcoor
   }

   draw_symmetry_elements
}

# Switch principle axes of inertia on/off
proc ::Symmetry::toggle_inertia {} {
   variable showinertia
   variable origmol
   variable idealmol
   variable com
   variable inertia
   if {![llength $inertia] || $origmol<0} { return }

   variable inertiagido
   variable inertiagidi
   foreach gid $inertiagido {
      graphics $origmol delete $gid
   }
   foreach gid $inertiagidi {
      graphics $idealmol delete $gid
   }

   if {$showinertia} {
      # Draw axes of inertia into the idealized coords molecule
      foreach a $inertia {
         lappend priaxes [lindex $a 0]
         lappend unique  [lindex $a 2]
      }
      puts "unique={$unique}"
      set all [atomselect $idealmol all]
      set inertiagidi [draw_axes_of_inertia $idealmol $all $com $priaxes $unique]

      # Draw axes of inertia into the original molecule
      variable sel
      foreach {rcom priaxes eigenval} [measure inertia $sel -eigenvals] {break}
      set eigenval [vecnorm $eigenval]

      set unique2 {0 0 0}
      if {[lindex $unique 0] && 
          abs([lindex $eigenval 1]-[lindex $eigenval 2])<0.03} {
         lset unique2 0 1
      }
      if {[lindex $unique 1] && 
          abs([lindex $eigenval 0]-[lindex $eigenval 2])<0.03} {
         lset unique2 1 1
      }
      if {[lindex $unique 2] && 
          abs([lindex $eigenval 0]-[lindex $eigenval 1])<0.03} {
         lset unique2 2 1
      }
      set inertiagido [draw_axes_of_inertia $origmol $sel $rcom $priaxes $unique2]
   }
}

# Determine the point group symmetry of a selection.
# The name of the point group is returned, possible values are
# C1, Ci, Cs, Cn, Cnv, Cinfv, Cnh, Sn, Dn, Dnh, Dnd, Dinfh,
# T, Td, Th, O, Oh, I, Ih, Kh, where n denotes the order of the highest 
# symmetry axis.
# Additionally you can get the symmetry elements (axes, mirror planes and 
# rotary-reflection axis) using the switch -elements.
# The algorithm is somewhat tolerant regarding perturbed coordinates,
# i.e. if the atoms don't sit exactly at the ideal position according
# to the symmetry operation, the point group will still be guessed.
# If the atoms are not too far off, the guess will be correct.
# The tolerance can be controlled through the -tol option.
#
# First the center of mass is computed through which all symmetry elements
# must pass. Diagonalizing the moments of inertia tensor we obtain the 
# moments of inertia as the eigenvalues and the principal axes of inertia
# as eigenvectors which contain some information about the symmetry.
# a) If one (only one) moment of inertia is zero and the other two have
#    the samevalue then the molecule must be linear. The point group is
#    Dinfh in case there is a center of inversion, Cinfn otherwise.
# b) If two of the moments of inertia are identical the molecule is a 
#    symmetric top. The possible points groups are C1, Cn, Cnv, Cnh, Sn,
#    Dn, Dnh, Dnd. The unique principal axis of inertia corresponds to one
#    of the rotation axes.
# c) If all three moments of inertia are different the molecule is an
#    asymmetric top. That means that the order of an axis cannot exceed 2.
#    Possible point groups are C1, Ci, Cs, C2, C2v, C2h, D2h, D2.
#    If at least one of the principal axes is C2 you have either D2h
#    (all principal axes are C2 and an inverson centre exists), D2 
#    (as D2h, but no inverson centre), C2v, C2h, or C2; otherwise you
#    have Ci, Cs, or C1.
# d) If all three moments of inertia are the same the molecule is a
#    spherical top. Possible point groups are T, Td, Th, O, Oh, I, Ih.
#    For spherical tops the rotation axes cannot be determined from the
#    principal axes because the latter are arbitrarily oriented.
#    However, every symmetry element has to pass through the center of
#    mass. Hence, you can calculate
#     the distance atom from the centre of gravity and if you find a set of atoms
#     of the same element with the same distances from the centre of gravity you
#     have found a rotation axis.
#     Note: There are a few cases of so-called accidiental spherical top molecules
#     which do not belong to a cubic point group. These are hard to handle, but
#     rare, too.

proc ::Symmetry::guess {sel {tol 0.2}} {
   variable origmol
   if {$origmol<0 || [catch {$sel num}]} { return }

#   variable checkbonds
#   if {$checkbonds} {
#      array set symm [measure symmetry $sel -tol $tol]
#   } else {
#      array set symm [measure symmetry $sel -tol $tol -nobonds]
#   }

   set allatom [atomselect top all]
   set m0 [transoffset [vecinvert [measure center $allatom weight mass]]]
   $allatom move $m0

   set vmdout [open "vmdout.xyz" w]
   set output [$sel get "type x y z"]
   set output [string map {"{" "" "}" "\n"} $output]
   puts $vmdout [$sel num]
   puts $vmdout " "
   puts $vmdout $output
   close $vmdout
   
   exec symm.exe vmdout.xyz tol=$tol
   
   #  "生成调用操作"
   set fp [open "a.tcl" r]
   set cmd [read  $fp]
   eval  $cmd ;#执行调用的文件
   close  $fp
   


   puts "pointgroup $pointgroup"
   puts "elements $elements"
   puts "com       $com"
   puts "inversion $inversion"
   puts "axes      $axes"
   puts "rraxes    $rraxes"
   puts "planes    $planes"
   
   puts [format "Pointgroup: %s, rmsd = %.2f" $pointgroup $rmsd]
   puts "Elements:   $elements"
#   if {[llength $symm(missing)]} {
#     puts "Missing:    $symm(missing)"
#   }
#   if {[llength $symm(additional)]} {
#     puts "Additional: $symm(additional)"
#   }

#   array unset symm

   variable uniquerepo
   if {[llength $uniquerepo]} {
      variable unique
      mol modselect [mol repindex $origmol $uniquerepo] $origmol "index $unique"
   }

   variable origcoor [$sel get {x y z}]
   variable showorient


   # Create molecule with idealized coordinates


   set filename [file rootname [molinfo $origmol get name]]
   
   save_viewpoint
   restore_viewpoint
   

   variable idealcoor

   variable showorient


   #toggle_ideal
   #toggle_inertia

   draw_symmetry_elements
   mol new vmdout.xyz
   mol top $origmol
}


proc ::Symmetry::draw_symmetry_elements {} {
   variable origmol
   variable sel
   if {$origmol<0 || [catch {$sel num}]} { return }

   # Setup empty molecules as drawing canvases
   save_viewpoint
   variable canvaso
   variable canvast
   if {$canvaso<0} {
      set canvaso [mol new]
      mol rename $canvaso "Canvas (opaque)"
   }
   if {$canvast<0} {
      set canvast [mol new]
      mol rename $canvast "Canvas (transparent)"
   }
   variable viewpoints
   set_viewpoint $viewpoints($origmol)

   graphics $canvaso delete all
   graphics $canvast delete all
   graphics $canvast material Transparent

   foreach {min max} [measure minmax $sel] {break}
   set radius [expr {0.5*[veclength [vecsub $min $max]]}]

   variable showorient
   variable orient
   variable com
   set center $com



   # Draw inversion center
   variable showcenter
   variable inversion
   variable centerinfo {}
   variable centerid {}

   if {$showcenter} {
       if {$inversion} {
           set ids [::Symmetry::draw_center $sel $center -label "  inv" -color purple]
           set centerid $ids
           lappend centerinfo "inv"
           variable lastselectcenter [lindex $centerid 0]
       }
   }


   variable showrotaxes
   variable axes
   if {$showrotaxes} {
      # Draw rotary axes
      variable axislist {}
      variable axisgid {}
      set i 0
      foreach axisobj $axes {
         foreach {axis order type} $axisobj {break}
         
         if {$order==-1} { set order inf }
         lappend axislist [format "%2i C%-3s %s" $i $order $type]
         
         if {$showorient} {
            set axis [coordtrans $rot $axis]
         }
         
         lappend axisgid [draw_symaxis $sel $center $axis -label "$i: C$order" \
                             -color orange3]
         incr i
      }
      variable lastselectaxis [lindex $axisgid 0 0]
   }

   variable showrraxes
   if {$showrraxes} {
      # Draw rotary reflection axes
      variable rraxes
      variable rraxislist {}
      variable rraxisgid {}
      set i 0
      foreach axisobj $rraxes {
         foreach {axis order type} $axisobj {break}

         if {$order==-1} { set order inf }

         lappend rraxislist [format "%2i S%-3s %s" $i $order $type]
         if {$showorient} {
            set axis [coordtrans $rot $axis]
         }

         lappend rraxisgid [draw_symaxis $sel $center $axis -label "S$order" \
                               -color lime]
         incr i
      }
      variable lastselectrraxis [lindex $rraxisgid 0 0]
   }

   variable showplanes
   if {$showplanes} {
      # Draw planes
      variable planes
      variable planelist {}
      variable planegid {}
      set i 0
      foreach planeobj $planes {
         foreach {normal type} $planeobj {break}
         lappend planelist [format "%2i %s" $i $type]
         
         if {$showorient} {
            set normal [coordtrans $rot $normal]
         }
         
         lappend planegid [draw_plane $center $normal -radius $radius \
                              -color pink]
         incr i
      }
      variable lastselectplaneo [lindex $planegid 0 0 0]
      variable lastselectplanet [lindex $planegid 0 1 0]
   }
   #puts "orient:"
   #mat_print $orient "% .4f"
}


proc ::Symmetry::save_viewpoint {} {
   variable viewpoints
   if [info exists viewpoints] {unset viewpoints}
   # get the current matricies
   foreach mol [molinfo list] {
      set viewpoints($mol) [molinfo $mol get {
        center_matrix rotate_matrix scale_matrix global_matrix}]
   }
   set top [molinfo top]
   return "$viewpoints($top)"
}

proc ::Symmetry::restore_viewpoint {} {
   variable viewpoints
   foreach mol [molinfo list] {
      if [info exists viewpoints($mol)] {
         molinfo $mol set {center_matrix rotate_matrix scale_matrix
           global_matrix} $viewpoints($mol)
      }
   }
}

proc ::Symmetry::set_viewpoint { viewpoint } {
   foreach mol [molinfo list] {
      molinfo $mol set {center_matrix rotate_matrix scale_matrix
           global_matrix} $viewpoint
   }
}

proc ::Symmetry::draw_symaxis {sel pivot axis args} {
   set color  [::util::getargs $args "-color" orange3]
   set label  [::util::getargs $args "-label" {}]
   set radius [::util::getargs $args "-radius" 0.05]

   set axis [vecnorm $axis]
   set minmax [measure minmax $sel]
   set diameter [veclength [vecsub [lindex $minmax 0] [lindex $minmax 1]]]

   set p1 [vecadd $pivot [vecscale $axis [expr {0.5*$diameter*1.1}]]]
   set p2 [vecsub $pivot [vecscale $axis [expr {0.5*$diameter*1.1}]]]
   set p3 [vecadd $pivot [vecscale $axis [expr {0.5*$diameter*1.15}]]]

   variable canvaso
   lappend ids [graphics $canvaso color $color]

   lappend ids [graphics $canvaso cylinder $p1 $p2 radius $radius resolution 11]
   lappend ids [graphics $canvaso sphere $p1 radius $radius resolution 11]
   lappend ids [graphics $canvaso sphere $p2 radius $radius resolution 11]
   if {[llength $label]} {
      lappend $ids [graphics $canvaso text $p3 $label]
   }

   return $ids
}

proc ::Symmetry::draw_plane {pivot normal args} {
   set radius   [::util::getargs $args "-radius" 1.0]
   set axis     [::util::getargs $args "-axis" {}]
   set nvplanes [::util::getargs $args "-nvplanes" 1]
   set color    [::util::getargs $args "-color" pink]
   set n [::util::getargs $args "-n" 50]

   set normal [vecnorm $normal]

   variable canvaso
   variable canvast
   lappend tids [graphics $canvast color $color]
   lappend oids [graphics $canvaso color $color]
   
   lappend tids [vmd_draw_disc $canvast $pivot $normal $radius \
		    n $n resolution 19]
   lappend oids [vmd_draw_ring $canvaso $pivot $normal $radius \
		    n $n width 0.02 resolution 19]

   for {set i 1} {$i<$nvplanes} {incr i} {
      set rot  [transabout $axis [expr {$i*360.0/$nvplanes}] deg]
      set plane [coordtrans $rot $normal]

      lappend tids [vmd_draw_disc $canvast $pivot $plane $radius n $n \
		       resolution 19]
      lappend oids [vmd_draw_ring $canvaso $pivot $plane $radius n $n \
		       width 0.02 resolution 19]
   }

   return [list $oids $tids]
}


proc ::Symmetry::draw_center {sel center args} {
    set color  [::util::getargs $args "-color" purple]
    set label  [::util::getargs $args "-label" {}]
    set radius [::util::getargs $args "-radius" 0.2]

    variable canvaso
    lappend ids [graphics $canvaso color $color]

    lappend ids [graphics $canvaso sphere $center radius $radius resolution 11]
    if {[llength $label]} {
        lappend ids [graphics $canvaso text $center $label]
    }

    return $ids
}

proc ::Symmetry::vmd_draw_priaxes { mol sel } {
   set COM [measure center $sel weight mass]
   
   foreach {com priaxes} [measure inertia $sel] {}

   draw_axes_of_inertia $mol $sel $COM $priaxes
}

proc ::Symmetry::draw_axes_of_inertia {mol sel COM priaxes {unique {}}} {
    # find the size of the system
    set minmax [measure minmax $sel]
    set ranges [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set scale [expr {0.95*[::util::lmax [lrange $ranges 0 2]]}]
    set scale2 [expr {1.02 * $scale}]
    set thin [expr {0.5*$scale}]

    foreach a $priaxes u $unique i {0 1 2} {
       # draw some nice vectors
       if {$u} {
          lappend ids [graphics $mol color white]
          lappend ids [draw_vector $mol $COM [vecscale $scale $a]]
       } else {
          lappend ids [graphics $mol color gray]
          lappend ids [draw_vector $mol $COM [vecscale $thin $a]]
       }
       lappend ids [graphics $mol text [vecadd $COM [vecscale $scale2 $a]] $i]
    }
    return [join $ids]
}

proc ::Symmetry::draw_scaled_arrow {mol start end args} {
   set res    [expr {   int([::util::getargs $args "resolution" 6])}]
   #set rad    [expr {double([::util::getargs $args "radius"  1.0])}]
   #set aspect [expr {double([::util::getargs $args "aspect" 5.0])}]

   set dir [vecsub $end $start]
   if {![veclength $dir]} { return }

   set scaling [expr [veclength $dir]/100]
   # an arrow is made of a cylinder and a cone
   set middle [vecadd $start [vecscale 0.8 $dir]]
   lappend ids [graphics $mol cylinder $start $middle radius [expr 2*$scaling] resolution $res]
   #puts [list cone $middle $end radius [expr 5*$scaling]]
   lappend ids [graphics $mol cone $middle $end radius [expr 5*$scaling] resolution $res]
   return $ids
}


proc ::Symmetry::draw_vector { mol pos val } {
   set end   [ vecadd $pos [ vecscale +1 $val ] ]
   return [draw_scaled_arrow $mol $pos $end]
}


# Draws a circle or a segement of a circle either using lines or tubes
# $center is the center of the circle, $bond is a vector perpendicular
# to the circle plane.
# By specifying ->, <-, or <-> aftre the keyword "arrow" you can add an
# arrow head to the end, the beginning or to both ends of the circle.
proc ::Symmetry::vmd_draw_ring {mol center bond ringrad args} {
   set resolution [expr {   int([::util::getargs $args "resolution" 6])}]
   set width      [expr {double([::util::getargs $args "width"  0.0])}]
   set offset     [expr {double([::util::getargs $args "offset"  0.0])}]
   set range      [expr {double([::util::getargs $args "range"  360.0])}]
   set aspect     [expr {double([::util::getargs $args "aspect" 5.0])}]
   set n          [expr {   int([::util::getargs $args "n"  10])}]
   set arrow      [::util::getargs $args "arrow" 0]
   set arrowsize  [expr {   int([::util::getargs $args "arrowsize" 6])}]
   set zero       [::util::getargs $args "zero" {0 1 0}]
   set end        [::util::getargs $args "end"  {}]
   set style      [::util::getargs $args "style" tube]
   set fill       [::util::getswitch $args "fill"]

   #graphics $mol color green
   #graphics $mol arrow $center [vecadd $center $zero]  radius $width resolution $resolution
   set coor1 [vecadd $center $bond]
   set cross [veccross $zero $bond]
   if {[veclength $cross]} {
      set zero [vecscale $ringrad [vecnorm $cross]]
   } else {
      # $zero and $bond are parallel
      lappend rnd [expr {rand()}]
      lappend rnd [expr {rand()}]
      lappend rnd [expr {rand()}]
      set zero [vecscale $ringrad [vecnorm [veccross $rnd $bond]]]
   }
   set zero [coordtrans [trans bond {0 0 0} $bond 90 deg] $zero]

   if {[llength $end]} {
      set range [vecangle $zero $end]
   }

   set arc $range
   set arrowhead 0.0
   set arrowoffset 0.0
   if {$arrow=="<-" || $arrow=="<->" || $arrow=="->"} {
      set arrowhead [expr {180.0/(3.1415*$ringrad)*($aspect)*$width}]
   }

   if {$arrow=="<-" || $arrow=="<->"} {
      set arc [expr {$arc-$arrowhead}]
      set arrowoffset $arrowhead
   }
   if {$arrow=="->" || $arrow=="<->"} {
      set arc [expr {$arc-$arrowhead}]
   }

   set m [trans bond $center $coor1 [expr {$offset+$arrowoffset}] deg]
   set start [coordtrans $m [vecadd $center $zero]]
   #graphics $mol color orange
   #graphics $mol arrow $center $start radius $width resolution $resolution
   #graphics $mol color red

   set oldp $start
   for {set i 1} {$i<=$n} {incr i} {
      set angle [expr {$arc*$i/$n}]
      set m [trans bond $center $coor1 $angle deg]
      set p [coordtrans $m $start]
      if {$style=="tube"} {
	 graphics $mol cylinder $oldp $p radius $width resolution $resolution
	 graphics $mol sphere $oldp radius $width resolution $resolution
      } else {
	 graphics $mol line $oldp $p width [expr {int($width)}]
      }
      if {$fill} {
	 graphics $mol triangle $center $oldp $p
      }
      set oldp $p
   }
   set m [trans bond $center $coor1 $arc deg]
   set p [coordtrans $m $start]
   if {$style=="tube"} {
      graphics $mol sphere $p radius $width resolution $resolution
   }

   if {$arrow=="->" || $arrow=="<->"} {
      set m [trans bond $center $coor1 [expr {$range-$arrowoffset}] deg]
      set tip [coordtrans $m $start]
      if {$style=="tube"} {
	 graphics $mol cone $p $tip radius [expr {2.2*$width}] resolution $resolution
      }
   } 
   if {$arrow=="<-" || $arrow=="<->"} {
      set m [trans bond $center $coor1 [expr {-$arrowhead}] deg]
      set tip [coordtrans $m $start]
      if {$style=="tube"} {
	 graphics $mol cone $start $tip radius [expr {2.2*$width}] resolution $resolution
      }
   }
}

proc ::Symmetry::vmd_draw_disc {mol center normal discrad args} {
   #set ringradius [expr {double([::util::getargs $args "ringradius"  0.0])}]

   eval vmd_draw_ring $mol [list $center] [list $normal] $discrad fill $args
}

