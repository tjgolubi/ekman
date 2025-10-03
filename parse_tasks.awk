BEGIN {
  field = "field0"
  fieldNum = 0
  innerNum = 0
  fname = "fname-" field "-" fieldNum ".wgs84"
}

/^ *<FRM A=/ {
  where = match($0, / A="([[:alnum:]]+)" B="([^"]+)"/, arr)
  if (where == 0) {
    print "Bad farm line: ", $0 > "/dev/stderr"
    next
  }
  name = arr[2]
  gsub(/ +/, "_", name)
  farms[arr[1]] = name
  next
}

/^ *<PFD A=".*" C="/ {
  where = match($0, / C="([^"]+)" D="[[:digit:]]+" E="[[:alnum:]]+" F="([[:alnum:]]+)">/, arr)
  if (where == 0) {
    print "Bad field line: ", $0 > "/dev/stderr"
    next
  }
  farm = farms[arr[2]]
  field = arr[1]
  gsub(/ +/, "_", field)
  fieldNum = 0
  innerNum = 0
  print farm field
  next
}

/^ *<LSG A="1".*>/ {
  ++fieldNum
  fname = farm "-" field "_" fieldNum ".wgs84"
  innerNum = 0
  print "" > fname
  print fname
  next
}

/^ *<LSG A="2".*>/ {
  print "# inner" > fname
  next
}

/^ *<LSG A="/ {
  print "Bad field section: ", $0 > "/dev/stderr"
  next
}

/^ *<PNT A="[[:digit:]]+" C="-?[0-9.]+" D="-?[0-9.]+"/ {
  where = match($0, / C="(-?[0-9.]+)" D="(-?[0-9.]+)"/, arr)
  if (where == 0) {
    print "Bad point line: ", $0 > "/dev/stderr"
    next
  }
  print arr[1], arr[2] >> fname
  next
}
