{{name}} <- structure(
  function() get("{{name}}"),
  class = "dust_system_generator",
  name = "{{name}}",
  properties = list(
    time_type = "{{time_type_property}}",
    has_compare = {{has_compare}},
    has_adjoint = {{has_adjoint}}),
  default_dt = {{default_dt}})
