within ;
package BCVTB_Graz_Test
  "Contains Model of 16 customer network connected to BCVTB system for cosimulation"

  package Components

    model Umrechnung
      "Control component which evaluates deviation from nominal mass flow rate"
      parameter Real m_flow_nominal;
      Modelica.Blocks.Interfaces.RealInput m_flow
      annotation (Placement(
            transformation(
            origin={10,106},
            extent={{-20,-20},{20,20}},
            rotation=270), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80})));
      Modelica.Blocks.Interfaces.RealOutput opening annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={2,-108}), iconTransformation(extent={{98,-10},{118,10}})));

    equation
      opening = m_flow/m_flow_nominal;
      annotation ();
    end Umrechnung;

    package Pipe

      package SpatialPipe

        model SpatialPipeVolume
          import BCVTB_Graz_Test;

          Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium
              = Medium)
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
          Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium
              = Medium)
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialMedium
            annotation (__Dymola_choicesAllMatching=true);
          parameter Modelica.SIunits.Length length=30 "Length of the pipe"
        annotation (Dialog(group="Geometry"));
          parameter Modelica.SIunits.Diameter diameter=0.1 "Pipe diameter"
        annotation (Dialog(group="Geometry"));
          parameter Modelica.SIunits.Length InsThickness=0.1
            "Thickness of the pipe insulation"
        annotation (Dialog(group="Geometry"));

          parameter Modelica.SIunits.ThermalConductivity lambda=0.026
            "Thermal conductivity of insulation"
              annotation (Dialog(group="Heat Loss"));
          parameter Modelica.SIunits.ThermalConductivity R=1/(lambda*2*
              Modelica.Constants.pi/Modelica.Math.log((diameter/2 + InsThickness)
              /(diameter/2)))
            "Thermal transmission coefficient between pipe medium and surrounding"
        annotation (Dialog(group="Heat Loss"));

          parameter Modelica.SIunits.MassFlowRate m_flow_nominal
            "Nominal mass flow rate"
          annotation(Dialog(group = "Nominal condition"));

          parameter Real V= ((diameter/2)^2)*Modelica.Constants.pi*length "Volume"
        annotation(Dialog(tab = "Advanced"));
          Modelica.SIunits.MassFlowRate m_flow "Mass flow in the pipe";

          parameter Modelica.SIunits.MassFlowRate m_flow_small(min=0) = 1E-4*abs(m_flow_nominal)
            "Small mass flow rate for regularization of zero flow"
          annotation(Dialog(tab = "Advanced"));
          parameter Modelica.SIunits.MassFlowRate m_flow_min = 0.0002
            "Minimum mass flow rate for regularization of zero flow"
          annotation(Dialog(tab = "Advanced"));

          BCVTB_Graz_Test.Components.Pipe.SpatialPipe.parts.flowResistance flowResistance(
            redeclare package Medium = Medium,
            diameter=diameter,
            length=length,
            redeclare function FlowModel = FlowModel)
            annotation (Placement(transformation(extent={{-48,-10},{-28,10}})));
          replaceable function FlowModel =
              Modelica.Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.massFlowRate_dp
              annotation(Dialog(tab = "Advanced"));

          BCVTB_Graz_Test.Components.Pipe.SpatialPipe.parts.spatialDelayInput heatTransportAndLoss(
            redeclare package Medium = Medium,
            length=length,
            T_amb=T_amb,
            lambda=lambda,
            R=R,
            InsThickness=InsThickness,
            diameter=diameter,
            m_flow=m_flow,
            T_start=T_start)
            annotation (Placement(transformation(extent={{0,-10},{20,10}})));
          BCVTB_Graz_Test.Components.Pipe.SpatialPipe.parts.ClosedVolume volume(
            nPorts=2,
            redeclare package Medium = Medium,
            use_portsData=false,
            V=V,
            scale=scale,
            p_start=p_start,
            use_T_start=use_T_start,
            T_start=T_start,
            h_start=h_start,
            X_start=X_start,
            C_start=C_start,
            use_HeatTransfer=false)
            annotation (Placement(transformation(extent={{62,2},{82,22}})));

          parameter Real scale=10;

          parameter Modelica.Media.Interfaces.Types.AbsolutePressure p_start=volume.system.p_start
            "Start value of pressure" annotation (Dialog(tab="VolumeInit"));
          parameter Boolean use_T_start=true "= true, use T_start, otherwise h_start"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.Temperature T_start=volume.system.T_start
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.SpecificEnthalpy h_start=if volume.use_T_start
               then Medium.specificEnthalpy_pTX(
              volume.p_start,
              volume.T_start,
              volume.X_start) else Medium.h_default "Start value of specific enthalpy"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.MassFraction X_start[Medium.nX]=
              Medium.X_default "Start value of mass fractions m_i/m"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.ExtraProperty C_start[Medium.nC]=
              fill(0, Medium.nC) "Start value of trace substances"
            annotation (Dialog(tab="VolumeInit"));

        public
          Modelica.Blocks.Interfaces.RealInput T_amb
            "Ambient temperature for pipe's surroundings" annotation (Placement(
                transformation(
                extent={{-20,-20},{20,20}},
                rotation=270,
                origin={0,102}), iconTransformation(
                extent={{-20,-20},{20,20}},
                rotation=270,
                origin={0,52})));
        equation
          m_flow= noEvent(if abs(port_a.m_flow) > m_flow_min then port_a.m_flow else m_flow_small);
          connect(port_a, flowResistance.port_a)
            annotation (Line(points={{-100,0},{-48,0}}, color={0,127,255}));
          connect(flowResistance.port_b, heatTransportAndLoss.port_a)
            annotation (Line(points={{-27.8,0},{0,0}},  color={0,127,255}));
          connect(heatTransportAndLoss.port_b, volume.ports[1])
            annotation (Line(points={{20,0},{46,0},{46,2},{70,2}},
                                                     color={0,127,255}));
          connect(port_b, volume.ports[2])
            annotation (Line(points={{100,0},{88,0},{88,2},{74,2}},
                                                      color={0,127,255}));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Rectangle(
                  extent={{-100,40},{100,-42}},
                  lineColor={0,0,0},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-100,32},{100,-34}},
                  lineColor={0,0,0},
                  fillColor={85,170,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-100,26},{100,-26}},
                  lineColor={0,0,0},
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={0,0,255})}),           Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
        end SpatialPipeVolume;

        model SpatialPipeVolume_DHN "Dual Flow Pipe Model"
          import DHNLibrary;

          Modelica.Fluid.Interfaces.FluidPort_a port_a1(redeclare package
              Medium =
                Medium)
            annotation (Placement(transformation(extent={{-110,30},{-90,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b port_a2(redeclare package
              Medium =
                Medium)
            annotation (Placement(transformation(extent={{90,-70},{110,-50}})));
          Modelica.Fluid.Interfaces.FluidPort_a port_b2(redeclare package
              Medium =
                Medium)
            annotation (Placement(transformation(extent={{-110,-70},{-90,-50}})));
          Modelica.Fluid.Interfaces.FluidPort_b port_b1(redeclare package
              Medium =
                Medium)
            annotation (Placement(transformation(extent={{90,30},{110,50}})));
          replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation (
              __Dymola_choicesAllMatching=true);

          parameter Modelica.SIunits.ThermalConductivity lambda=0.026
            "Thermal conductivity of insulation"
            annotation (Dialog(group="Heat Loss"));

          parameter Modelica.SIunits.Length length=30 "Length of the pipe"
        annotation (Dialog(group="Geometry"));
          parameter Modelica.SIunits.Diameter diameter=0.1 "Pipe diameter"
        annotation (Dialog(group="Geometry"));
          parameter Modelica.SIunits.Length InsThickness=0.1
            "Thickness of the pipe insulation"
        annotation (Dialog(group="Geometry"));
          parameter Modelica.SIunits.MassFlowRate m_flow_nominal=10
            "Nominal mass flow rate"
            annotation (Dialog(group="Nominal Values"));
          SpatialPipeVolume                                                            pipe(
            length=length,
            diameter=diameter,
            InsThickness=InsThickness,
            lambda=lambda,
            R=R,
            m_flow_nominal=m_flow_nominal,
            V=V,
            redeclare package Medium = Medium,
            scale=scale,
            m_flow_small=m_flow_small,
            m_flow_min=m_flow_min,
            redeclare function FlowModel = FlowModel,
            p_start=p_start,
            use_T_start=use_T_start,
            T_start=T_start,
            h_start=h_start,
            X_start=X_start,
            C_start=C_start)
            annotation (Placement(transformation(extent={{-26,22},{26,58}})));
          SpatialPipeVolume                                                            pipe1(
            length=length,
            diameter=diameter,
            InsThickness=InsThickness,
            lambda=lambda,
            R=R,
            m_flow_nominal=m_flow_nominal,
            V=V,
            redeclare package Medium = Medium,
            scale=scale,
            m_flow_small=m_flow_small,
            m_flow_min=m_flow_min,
            redeclare function FlowModel = FlowModel,
            p_start=p_start,
            use_T_start=use_T_start,
            T_start=T_start,
            h_start=h_start,
            X_start=X_start,
            C_start=C_start)
            annotation (Placement(transformation(extent={{-26,-80},{24,-40}})));
          parameter Modelica.Media.Interfaces.Types.AbsolutePressure p_start=
              pipe.volume.system.p_start "Start value of pressure"
         annotation (Dialog(tab="VolumeInit"));
          parameter Real scale=10;
          parameter Modelica.SIunits.ThermalConductivity R=1/(pipe.lambda*2*Modelica.Constants.pi
              /Modelica.Math.log((pipe.diameter/2 + pipe.InsThickness)/(pipe.diameter/2)))
            "Thermal transmission coefficient between pipe medium and surrounding"
        annotation (Dialog(group="Heat Loss"));

          parameter Boolean use_T_start=true "= true, use T_start, otherwise h_start"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.Temperature T_start=pipe.volume.system.T_start
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.SpecificEnthalpy h_start=if pipe.volume.use_T_start
               then Medium.specificEnthalpy_pTX(
              pipe.volume.p_start,
              pipe.volume.T_start,
              pipe.volume.X_start) else Medium.h_default
            "Start value of specific enthalpy" annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.MassFraction X_start[Medium.nX]=
              Medium.X_default "Start value of mass fractions m_i/m"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.Media.Interfaces.Types.ExtraProperty C_start[Medium.nC]=
              fill(0, Medium.nC) "Start value of trace substances"
            annotation (Dialog(tab="VolumeInit"));
          parameter Modelica.SIunits.MassFlowRate m_flow_small=1E-4*abs(pipe.m_flow_nominal)
            "Small mass flow rate for regularization of zero flow"
            annotation (Dialog(tab="Advanced"));
          parameter Modelica.SIunits.MassFlowRate m_flow_min=0.0002
            "Minimum mass flow rate for regularization of zero flow"
            annotation (Dialog(tab="Advanced"));
          replaceable function FlowModel =
              Modelica.Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.massFlowRate_dp
            annotation (Dialog(tab="Advanced"));
          parameter Modelica.SIunits.SpecificEnthalpy initialEnthalpy[:]={Medium.h_default,
              Medium.h_default}
            "Initial values of the Enthalpy for the spatialDistribuiton"
            annotation (Dialog(tab="Advanced"));
          parameter Real V=((pipe.diameter/2)^2)*Modelica.Constants.pi*pipe.length "Volume"
                     annotation (Dialog(tab="Advanced"));
          Modelica.Blocks.Sources.RealExpression realExpression(y=260)
            annotation (Placement(transformation(extent={{-118,108},{-98,128}})));
        equation
          connect(port_a1, pipe.port_a)
            annotation (Line(points={{-100,40},{-26,40}}, color={0,127,255}));
          connect(pipe.port_b, port_b1) annotation (Line(points={{26,40},{26,40},{
                  100,40}}, color={0,127,255}));
          connect(pipe1.port_a, port_b2)
            annotation (Line(points={{-26,-60},{-100,-60}}, color={0,127,255}));
          connect(pipe1.port_b, port_a2) annotation (Line(points={{24,-60},{62,-60},
                  {100,-60}}, color={0,127,255}));
          connect(realExpression.y, pipe.T_amb) annotation (Line(points={{-97,118},
                  {-46,118},{-46,49.36},{0,49.36}}, color={0,0,127}));
          connect(realExpression.y, pipe1.T_amb) annotation (Line(points={{-97,118},
                  {-48,118},{-48,-49.6},{-1,-49.6}}, color={0,0,127}));
          annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Rectangle(
                  extent={{-100,80},{100,-2}},
                  lineColor={0,0,0},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-100,72},{100,6}},
                  lineColor={0,0,0},
                  fillColor={85,170,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-102,66},{98,14}},
                  lineColor={0,0,0},
                  fillColor={35,188,18},
                  fillPattern=FillPattern.CrossDiag),
                          Rectangle(
                  extent={{-100,-14},{100,-96}},
                  lineColor={0,0,0},
                  fillColor={0,0,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-100,-22},{100,-88}},
                  lineColor={0,0,0},
                  fillColor={85,170,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-100,-28},{100,-80}},
                  lineColor={0,0,0},
                  fillColor={35,188,18},
                  fillPattern=FillPattern.CrossDiag)}),                  Diagram(coordinateSystem(
                  preserveAspectRatio=false)));
        end SpatialPipeVolume_DHN;

        package parts

          connector FluidPorts_b
            "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
            extends Modelica.Fluid.Interfaces.FluidPort;
             // parameter Integer nPorts=0 "Number of ports" annotation(Dialog(connectorSizing=true));
            annotation (defaultComponentName="ports_b",
                        Diagram(coordinateSystem(
                  preserveAspectRatio=false,
                  extent={{-50,-200},{50,200}},
                  initialScale=0.2), graphics={
                  Text(extent={{-75,130},{75,100}}, textString="%name"),
                  Rectangle(
                    extent={{-25,100},{25,-100}},
                    lineColor={0,0,0}),
                  Ellipse(
                    extent={{-25,90},{25,40}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-25,25},{25,-25}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-25,-40},{25,-90}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-15,-50},{15,-80}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-15,15},{15,-15}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-15,50},{15,80}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid)}),
                 Icon(coordinateSystem(
                  preserveAspectRatio=false,
                  extent={{-50,-200},{50,200}},
                  initialScale=0.2), graphics={
                  Rectangle(
                    extent={{-50,200},{50,-200}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-50,180},{50,80}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-50,50},{50,-50}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-50,-80},{50,-180}},
                    lineColor={0,0,0},
                    fillColor={0,127,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-30,30},{30,-30}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-30,100},{30,160}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid),
                  Ellipse(
                    extent={{-30,-100},{30,-160}},
                    lineColor={0,127,255},
                    fillColor={255,255,255},
                    fillPattern=FillPattern.Solid)}));
          end FluidPorts_b;

          model ClosedVolume
            "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
          import Modelica.Constants.pi;

            // Mass and energy balance, ports
            extends PartialLumpedVessel(
              final fluidVolume=V,
              vesselArea=pi*(3/4*V)^(2/3),
              heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));

            parameter Modelica.SIunits.Volume V "Volume";

          equation
            Wb_flow = 0;
            for i in 1:nPorts loop
              vessel_ps_static[i] = medium.p;
            end for;

            annotation (defaultComponentName="volume",
              Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
                    100,100}}), graphics={Ellipse(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,0},
                  fillPattern=FillPattern.Sphere,
                  fillColor={170,213,255}), Text(
                  extent={{-150,12},{150,-18}},
                  lineColor={0,0,0},
                  textString="V=%V")}),
            Documentation(info="<html>
<p>
Ideally mixed volume of constant size with two fluid ports and one medium model.
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>.
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected.
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</p>
<p>
If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between volume and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>.
</p>
</html>"));
          end ClosedVolume;

          partial model PartialLumpedVessel
            "Lumped volume with a vector of fluid ports and replaceable heat transfer model"
            extends PartialLumpedVolume;

            // Port definitions
            parameter Integer nPorts=0 "Number of ports"
              annotation(Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

            Modelica.Fluid.Vessels.BaseClasses.VesselFluidPorts_b ports[nPorts](
                redeclare each package Medium = Medium) "Fluid inlets and outlets"
              annotation (Placement(transformation(extent={{-40,-10},{40,10}}, origin={0,-100})));

            // Port properties
            parameter Boolean use_portsData=true
              "= false to neglect pressure loss and kinetic energy"
              annotation(Evaluate=true, Dialog(tab="General",group="Ports"));
            parameter Modelica.Fluid.Vessels.BaseClasses.VesselPortsData[if use_portsData then nPorts else 0]
            portsData "Data of inlet/outlet ports"
              annotation(Dialog(tab="General",group="Ports",enable= use_portsData));

            parameter Medium.MassFlowRate m_flow_nominal = if system.use_eps_Re then system.m_flow_nominal else 1e2*system.m_flow_small
              "Nominal value for mass flow rates in ports"
              annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));
            parameter Modelica.SIunits.MassFlowRate m_flow_small(min=0)=if system.use_eps_Re
               then system.eps_m_flow*m_flow_nominal else system.m_flow_small
              "Regularization range at zero mass flow rate" annotation (Dialog(
                tab="Advanced",
                group="Port properties",
                enable=stiffCharacteristicForEmptyPort));
            parameter Boolean use_Re = system.use_eps_Re
              "= true, if turbulent region is defined by Re, otherwise by m_flow_small"
              annotation(Dialog(tab="Advanced", group="Port properties"), Evaluate=true);
          /*
  parameter Medium.AbsolutePressure dp_small = system.dp_small
    "Turbulent flow if |dp| >= dp_small (regularization of zero flow)"
    annotation(Dialog(tab="Advanced",group="Ports"));
*/
            Medium.EnthalpyFlowRate ports_H_flow[nPorts];
            Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
            Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
              "Substance mass flows through ports";
            Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];
            Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
              "Trace substance mass flows through ports";

            // Heat transfer through boundary
            parameter Boolean use_HeatTransfer = false
              "= true to use the HeatTransfer model"
                annotation (Dialog(tab="Assumptions", group="Heat transfer"));
            replaceable model HeatTransfer =
              Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
              constrainedby
              Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
              "Wall heat transfer"
                annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);

            HeatTransfer heatTransfer(
              redeclare final package Medium = Medium,
              final n=1,
              final states = {medium.state},
              final use_k = use_HeatTransfer)
                annotation (Placement(transformation(
                  extent={{-10,-10},{30,30}},
                  rotation=90,
                  origin={-50,-10})));
            Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

            // Conservation of kinetic energy
            Medium.Density[nPorts] portInDensities
              "densities of the fluid at the device boundary";
            Modelica.SIunits.Velocity[nPorts] portVelocities
              "velocities of fluid flow at device boundary";
            Modelica.SIunits.EnergyFlowRate[nPorts] ports_E_flow
              "flow of kinetic and potential energy at device boundary";

            // Note: should use fluidLevel_start - portsData.height
            Real[nPorts] s(each start = fluidLevel_max)
              "curve parameters for port flows vs. port pressures; for further details see, Modelica Tutorial: Ideal switching devices";
            Real[nPorts] ports_penetration
              "penetration of port with fluid, depending on fluid level and port diameter";

            // treatment of pressure losses at ports
            Modelica.SIunits.Area[nPorts] portAreas={Modelica.Constants.pi/4*
                portsData_diameter[i]^2 for i in 1:nPorts};
            Medium.AbsolutePressure[nPorts] vessel_ps_static
              "static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

            // determination of turbulent region
            constant Modelica.SIunits.ReynoldsNumber Re_turbulent=100
              "cf. suddenExpansion";
            Modelica.SIunits.MassFlowRate[nPorts] m_flow_turbulent;

          protected
            input Modelica.SIunits.Height fluidLevel=0
              "level of fluid in the vessel for treating heights of ports";
            parameter Modelica.SIunits.Height fluidLevel_max=1
              "maximum level of fluid in the vessel";
            parameter Modelica.SIunits.Area vesselArea=Modelica.Constants.inf
              "Area of the vessel used to relate to cross flow area of ports";

            // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
            // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
            // providing portsData_diameter and portsData_height, independent of the use_portsData setting.
            // Note: this moreover serves as work-around if a tool does not support a zero sized portsData record.
            Modelica.Blocks.Interfaces.RealInput[nPorts]
            portsData_diameter_internal = portsData.diameter if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal = portsData.height if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal = portsData.zeta_in if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out_internal = portsData.zeta_out if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;
            Modelica.Blocks.Interfaces.BooleanInput[nPorts] regularFlow(each start=true);
            Modelica.Blocks.Interfaces.BooleanInput[nPorts] inFlow(each start=false);

          equation
            mb_flow = sum(ports.m_flow);
            mbXi_flow = sum_ports_mXi_flow;
            mbC_flow  = sum_ports_mC_flow;
            Hb_flow = sum(ports_H_flow) + sum(ports_E_flow);
            Qb_flow = heatTransfer.Q_flows[1];

            // Only one connection allowed to a port to avoid unwanted ideal mixing
            for i in 1:nPorts loop
              assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");         end for;
            // Check for correct solution
            assert(fluidLevel <= fluidLevel_max, "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(fluidLevel) + ")");
            assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(fluidLevel) + ") is below zero meaning that the solution failed.");

            // Boundary conditions

            // treatment of conditional portsData
            connect(portsData_diameter, portsData_diameter_internal);
            connect(portsData_height, portsData_height_internal);
            connect(portsData_zeta_in, portsData_zeta_in_internal);
            connect(portsData_zeta_out, portsData_zeta_out_internal);
            if not use_portsData then
              portsData_diameter = zeros(nPorts);
              portsData_height = zeros(nPorts);
              portsData_zeta_in = zeros(nPorts);
              portsData_zeta_out = zeros(nPorts);
            end if;

            // actual definition of port variables
            for i in 1:nPorts loop
              portInDensities[i] = Medium.density(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)));
              if use_portsData then
                // dp = 0.5*zeta*d*v*|v|
                // Note: assume vessel_ps_static for portVelocities to avoid algebraic loops for ports.p
                portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/Medium.density(Medium.setState_phX(vessel_ps_static[i], actualStream(ports[i].h_outflow), actualStream(ports[i].Xi_outflow))));
                // Note: the penetration should not go too close to zero as this would prevent a vessel from running empty
                ports_penetration[i] = Modelica.Fluid.Utilities.regStep(
                  fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i],
                  1,
                  1e-3,
                  0.1*portsData_diameter[i]);
                m_flow_turbulent[i]=if not use_Re then m_flow_small else
                  max(m_flow_small, (Modelica.Constants.pi/8)*portsData_diameter[i]
                                     *(Medium.dynamicViscosity(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)))
                                       + Medium.dynamicViscosity(medium.state))*Re_turbulent);
              else
                // an infinite port diameter is assumed
                portVelocities[i] = 0;
                ports_penetration[i] = 1;
                m_flow_turbulent[i] = Modelica.Constants.inf;
              end if;

              // fluid flow through ports
              regularFlow[i] = fluidLevel >= portsData_height[i];
              inFlow[i]      = not regularFlow[i] and (s[i] > 0 or portsData_height[i] >= fluidLevel_max);
              if regularFlow[i] then
                // regular operation: fluidLevel is above ports[i]
                // Note: >= covers default values of zero as well
                if use_portsData then
                  /* Without regularization
                 ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2
                              * noEvent(if ports[i].m_flow>0 then zeta_in[i]/portInDensities[i] else -zeta_out[i]/medium.d);
              */

                  ports[i].p = homotopy(vessel_ps_static[i] + (0.5/portAreas[i]^2*
                    Modelica.Fluid.Utilities.regSquare2(
                    ports[i].m_flow,
                    m_flow_turbulent[i],
                    (portsData_zeta_in[i] - 1 + portAreas[i]^2/vesselArea^2)/
                      portInDensities[i]*ports_penetration[i],
                    (portsData_zeta_out[i] + 1 - portAreas[i]^2/vesselArea^2)/medium.d/
                      ports_penetration[i])), vessel_ps_static[i]);
                  /*
                // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
                ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                       2*portInDensities[i]/portsData_zeta_in[i],
                                       2*medium.d/portsData_zeta_out[i]));
              */
                else
                  ports[i].p = vessel_ps_static[i];
                end if;
                s[i] = fluidLevel - portsData_height[i];

              elseif inFlow[i] then
                // ports[i] is above fluidLevel and has inflow
                ports[i].p = vessel_ps_static[i];
                s[i] = ports[i].m_flow;

              else
                // ports[i] is above fluidLevel, preventing outflow
                ports[i].m_flow = 0;
                s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(portsData_height[i] - fluidLevel);
              end if;

              ports[i].h_outflow  = medium.h;
              ports[i].Xi_outflow = medium.Xi;
              ports[i].C_outflow  = C;

              ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
                "Enthalpy flow";
              ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[i] + system.g*portsData_height[i])
                "Flow of kinetic and potential energy";
              ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
                "Component mass flow";
              ports_mC_flow[i,:]  = ports[i].m_flow * actualStream(ports[i].C_outflow)
                "Trace substance mass flow";
            end for;

            for i in 1:Medium.nXi loop
              sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
            end for;

            for i in 1:Medium.nC loop
              sum_ports_mC_flow[i]  = sum(ports_mC_flow[:,i]);
            end for;

            connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
                points={{-100,0},{-87,0},{-87,8.88178e-016},{-74,8.88178e-016}},
                color={191,0,0},
                smooth=Smooth.None));
           annotation (
            Documentation(info="<html>
<p>
This base class extends PartialLumpedVolume with a vector of fluid ports and a replaceable wall HeatTransfer model.
</p>
<p>
The following modeling assumption are made:
<ul>
<li>homogeneous medium, i.e., phase separation is not taken into account,</li>
<li>no kinetic energy in the fluid, i.e., kinetic energy dissipates into the internal energy,</li>
<li>pressure loss definitions at vessel ports assume incompressible fluid,</li>
<li>outflow of ambient media is prevented at each port assuming check valve behavior.
    If <code> fluidlevel &lt; portsData_height[i] </code>and &nbsp; <code> ports[i].p &lt; vessel_ps_static[i]</code> massflow at the port is set to 0.</li>
</ul>
<p>
Each port has a (hydraulic) diameter and a height above the bottom of the vessel, which can be configured using the &nbsp;<b><code>portsData</code></b> record.
Alternatively the impact of port geometries can be neglected with <code>use_portsData=false</code>. This might be useful for early
design studies. Note that this means to assume an infinite port diameter at the bottom of the vessel.
Pressure drops and heights of the ports as well as kinetic and potential energy fluid entering or leaving the vessel are neglected then.
</p>
<p>
The following variables need to be defined by an extending model:
</p>
<ul>
<li><code>input fluidVolume</code>, the volume of the fluid in the vessel,</li>
<li><code>vessel_ps_static[nPorts]</code>, the static pressures inside the vessel at the height of the corresponding ports, at zero flow velocity, and</li>
<li><code>Wb_flow</code>, work term of the energy balance, e.g., p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
<p>
An extending model should define:
</p>
<ul>
<li><code>parameter vesselArea</code> (default: Modelica.Constants.inf m2), the area of the vessel, to be related to cross flow areas of the ports for the consideration of dynamic pressure effects.</li>
</ul>
<p>
Optionally the fluid level may vary in the vessel, which effects the flow through the ports at configurable <code>portsData_height[nPorts]</code>.
This is why an extending model with varying fluid level needs to define:
</p>
<ul>
<li><code>input fluidLevel (default: 0m)</code>, the level the fluid in the vessel, and</li>
<li><code>parameter fluidLevel_max (default: 1m)</code>, the maximum level that must not be exceeded. Ports at or above fluidLevel_max can only receive inflow.</li>
</ul>
<p>
An extending model should not access the <code>portsData</code> record defined in the configuration dialog,
as an access to <code>portsData</code> may fail for <code>use_portsData=false</code> or <code>nPorts=0</code>.
</p>
<p>
Instead the predefined variables
</p>
<ul>
<li><code>portsData_diameter[nPorts]</code>,</li>
<li><code>portsData_height[nPorts]</code>,</li>
<li><code>portsData_zeta_in[nPorts]</code>, and</li>
<li><code>portsData_zeta_out[nPorts]</code></li>
</ul>
<p>
should be used if these values are needed.
</p>
</html>",           revisions="<html>
<ul>
<li><i>Jan. 2009</i> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic and potential energy of fluid entering or leaving in energy balance</li>
   </ul>
</li>
<li><i>Dec. 2008</i> by R&uuml;diger Franke: derived from OpenTank, in order to make general use of configurable port diameters</li>
</ul>
</html>"),    Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                    {100,100}}), graphics={Text(
                  extent={{-150,110},{150,150}},
                  textString="%name",
                  lineColor={0,0,255})}));
          end PartialLumpedVessel;

          partial model PartialLumpedVolume "Lumped volume with mass and energy balance"
          import Modelica.Fluid.Types;
          import Modelica.Fluid.Types.Dynamics;
          import Modelica.Media.Interfaces.Choices.IndependentVariables;

            outer Modelica.Fluid.System system "System properties";
            replaceable package Medium =
              Modelica.Media.Interfaces.PartialMedium "Medium in the component"
                annotation (choicesAllMatching = true);

            // Inputs provided to the volume model
            input Modelica.SIunits.Volume fluidVolume "Volume";
            // scaling parameter for m_flow
            parameter Real scale = 10;
            // Assumptions
            parameter Types.Dynamics energyDynamics=system.energyDynamics
              "Formulation of energy balance"
              annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
            parameter Types.Dynamics massDynamics=system.massDynamics
              "Formulation of mass balance"
              annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
            final parameter Types.Dynamics substanceDynamics=massDynamics
              "Formulation of substance balance"
              annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
            final parameter Types.Dynamics traceDynamics=massDynamics
              "Formulation of trace substance balance"
              annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

            // Initialization
            parameter Medium.AbsolutePressure p_start = system.p_start
              "Start value of pressure"
              annotation(Dialog(tab = "Initialization"));
            parameter Boolean use_T_start = true "= true, use T_start, otherwise h_start"
              annotation(Dialog(tab = "Initialization"), Evaluate=true);
            parameter Medium.Temperature T_start=system.T_start
              annotation(Dialog(tab = "Initialization", enable = use_T_start));
                                                               /*
    if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start) 
    "Start value of temperature"*/
            parameter Medium.SpecificEnthalpy h_start=
              if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default
              "Start value of specific enthalpy"
              annotation(Dialog(tab = "Initialization", enable = not use_T_start));
            parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
              "Start value of mass fractions m_i/m"
              annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
            parameter Medium.ExtraProperty C_start[Medium.nC](
                 quantity=Medium.extraPropertiesNames)=fill(0, Medium.nC)
              "Start value of trace substances"
              annotation (Dialog(tab="Initialization", enable=Medium.nC > 0));

            Medium.BaseProperties medium(
              preferredMediumStates=true,
              p(start=p_start),
              h(start=h_start),
              T(start=T_start),
              Xi(start=X_start[1:Medium.nXi]));
            Modelica.SIunits.Energy U "Internal energy of fluid";
            Modelica.SIunits.Mass m "Mass of fluid";
            Modelica.SIunits.Mass[Medium.nXi] mXi
              "Masses of independent components in the fluid";
            Modelica.SIunits.Mass[Medium.nC] mC "Masses of trace substances in the fluid";
            // C need to be added here because unlike for Xi, which has medium.Xi,
            // there is no variable medium.C
            Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";

            // variables that need to be defined by an extending class
            Modelica.SIunits.MassFlowRate mb_flow "Mass flows across boundaries";
            Modelica.SIunits.MassFlowRate[Medium.nXi] mbXi_flow
              "Substance mass flows across boundaries";
            Medium.ExtraPropertyFlowRate[Medium.nC] mbC_flow
              "Trace substance mass flows across boundaries";
            Modelica.SIunits.EnthalpyFlowRate Hb_flow
              "Enthalpy flow across boundaries or energy source/sink";
            Modelica.SIunits.HeatFlowRate Qb_flow
              "Heat flow across boundaries or energy source/sink";
            Modelica.SIunits.Power Wb_flow "Work flow across boundaries or source term";
          protected
            parameter Boolean initialize_p = not Medium.singleState
              "= true to set up initial equations for pressure";
            Real[Medium.nC] mC_scaled(min=fill(Modelica.Constants.eps, Medium.nC))
              "Scaled masses of trace substances in the fluid";
          equation
            assert(not (energyDynamics<>Dynamics.SteadyState and massDynamics==Dynamics.SteadyState) or Medium.singleState,
                   "Bad combination of dynamics options and Medium not conserving mass if fluidVolume is fixed.");

            // Total quantities
            m = fluidVolume*medium.d;
            mXi = m*medium.Xi;
            U = m*medium.u;
            mC = m*C;

            // Energy and mass balances
            if energyDynamics == Dynamics.SteadyState then
              0 = Hb_flow + Qb_flow + Wb_flow;
            else
              der(U) = Hb_flow + Qb_flow + Wb_flow;
            end if;

            if massDynamics == Dynamics.SteadyState then
              0 = mb_flow;
            else
              der(m) = scale * mb_flow;
            end if;

            if substanceDynamics == Dynamics.SteadyState then
              zeros(Medium.nXi) = mbXi_flow;
            else
              der(mXi) = mbXi_flow;
            end if;

            if traceDynamics == Dynamics.SteadyState then
              zeros(Medium.nC)  = mbC_flow;
            else
              der(mC_scaled) = mbC_flow./Medium.C_nominal;
            end if;
              mC = mC_scaled.*Medium.C_nominal;

          initial equation
            // initialization of balances
            if energyDynamics == Dynamics.FixedInitial then
              /*
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    */
              if Medium.ThermoStates == IndependentVariables.ph or
                 Medium.ThermoStates == IndependentVariables.phX then
                 medium.h = h_start;
              else
                 medium.T = T_start;
              end if;
            elseif energyDynamics == Dynamics.SteadyStateInitial then
              /*
    if use_T_start then
      der(medium.T) = 0;
    else
      der(medium.h) = 0;
    end if;
    */
              if Medium.ThermoStates == IndependentVariables.ph or
                 Medium.ThermoStates == IndependentVariables.phX then
                 der(medium.h) = 0;
              else
                 der(medium.T) = 0;
              end if;
            end if;

            if massDynamics == Dynamics.FixedInitial then
              if initialize_p then
                medium.p = p_start;
              end if;
            elseif massDynamics == Dynamics.SteadyStateInitial then
              if initialize_p then
                der(medium.p) = 0;
              end if;
            end if;

            if substanceDynamics == Dynamics.FixedInitial then
              medium.Xi = X_start[1:Medium.nXi];
            elseif substanceDynamics == Dynamics.SteadyStateInitial then
              der(medium.Xi) = zeros(Medium.nXi);
            end if;

            if traceDynamics == Dynamics.FixedInitial then
              mC_scaled = m*C_start[1:Medium.nC]./Medium.C_nominal;
            elseif traceDynamics == Dynamics.SteadyStateInitial then
              der(mC_scaled) = zeros(Medium.nC);
            end if;

            annotation (
              Documentation(info="<html>
<p>
Interface and base class for an ideally mixed fluid volume with the ability to store mass and energy.
The following boundary flow and source terms are part of the energy balance and must be specified in an extending class:
</p>
<ul>
<li><code><b>Qb_flow</b></code>, e.g., convective or latent heat flow rate across segment boundary, and</li>
<li><code><b>Wb_flow</b></code>, work term, e.g., p*der(fluidVolume) if the volume is not constant.</li>
</ul>
<p>
The component volume <code><b>fluidVolume</b></code> is an input that needs to be set in the extending class to complete the model.
</p>
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><code><b>Hb_flow</b></code>, enthalpy flow,</li>
<li><code><b>mb_flow</b></code>, mass flow,</li>
<li><code><b>mbXi_flow</b></code>, substance mass flow, and</li>
<li><code><b>mbC_flow</b></code>, trace substance mass flow.</li>
</ul>
</html>"));
          end PartialLumpedVolume;

          model flowResistance

            Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package
                Medium =                                                            Medium)
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
            Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package
                Medium =                                                            Medium)
              annotation (Placement(transformation(extent={{92,-10},{112,10}})));

           replaceable function FlowModel =
                Modelica.Fluid.Pipes.BaseClasses.WallFriction.QuadraticTurbulent.massFlowRate_dp;

            replaceable package Medium =
                Modelica.Media.Interfaces.PartialMedium annotation (
                __Dymola_choicesAllMatching=true);

                parameter Modelica.SIunits.Diameter diameter;
                parameter Modelica.SIunits.Length length;
                Modelica.SIunits.Density rho;
                Modelica.SIunits.DynamicViscosity mu;
                Modelica.SIunits.Pressure dp;

          equation
            dp = port_a.p - port_b.p;
            mu = Medium.dynamicViscosity(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
            rho = Medium.density(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));

            port_a.m_flow = FlowModel(dp, rho, rho, mu, mu, length, diameter);

            0 = port_a.m_flow + port_b.m_flow;

            port_a.Xi_outflow = inStream(port_b.Xi_outflow);
            port_b.Xi_outflow = inStream(port_a.Xi_outflow);
            port_a.C_outflow = inStream(port_b.C_outflow);
            port_b.C_outflow = inStream(port_a.C_outflow);

            port_b.h_outflow = inStream(port_a.h_outflow);
            port_a.h_outflow = inStream(port_b.h_outflow);
            annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
                    extent={{-100,40},{100,-40}},
                    lineColor={28,108,200},
                    fillColor={0,0,255},
                    fillPattern=FillPattern.Forward), Line(
                    points={{-30,-28},{-30,32},{18,-28},{18,32}},
                    color={0,0,0},
                    thickness=1)}),                                        Diagram(
                  coordinateSystem(preserveAspectRatio=false)));
          end flowResistance;

          model spatialDelayInput

            Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package
                Medium =
                  Medium)
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
            Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package
                Medium =
                  Medium)
              annotation (Placement(transformation(extent={{90,-10},{110,10}})));
          replaceable package Medium =
                Modelica.Media.Interfaces.PartialMedium annotation (
                __Dymola_choicesAllMatching=true);
            constant Real pi=Modelica.Constants.pi;

            parameter Modelica.SIunits.Length length = 30 "Length of the pipe"
            annotation (Dialog(group="Geometry"));

            Modelica.SIunits.Volume Volume "Volume of the pipe";
            Modelica.SIunits.Density density;
            Modelica.SIunits.Time delayTime;
            Modelica.SIunits.Velocity velocity "Flow velocity of medium in pipe";

            //Used for second heatloss spatial
            Modelica.SIunits.SpecificEnthalpy h_a;
            Modelica.SIunits.SpecificEnthalpy h_b;
            Modelica.SIunits.Time tau;
            Modelica.SIunits.Time tauExp;
            Modelica.SIunits.Time time_b;
            Modelica.SIunits.Time time_a;

            Real exponential(min=0) "exponential-decay function";
            Modelica.SIunits.Temp_K T_in "Temperature of the incoming water";
            Modelica.SIunits.Temp_K T_out "Temperature of the outgoing water";
            Modelica.SIunits.SpecificHeatCapacity cp;
            input Modelica.SIunits.Temp_K T_amb = 273.15;
            parameter Modelica.SIunits.ThermalConductivity lambda=0.026
              "Thermal conductivity of insulation";
            parameter Modelica.SIunits.ThermalConductivity R = 1 / (lambda*2*Modelica.Constants.pi/Modelica.Math.log((diameter/2+InsThickness)/(diameter/2)))
              "Thermal transmission coefficient between pipe medium and surrounding";
            parameter Modelica.SIunits.Length InsThickness = 0.1
              "Thickness of the pipe insulation";

            parameter Modelica.SIunits.Diameter diameter = 0.1 "Pipe diameter"
          annotation (Dialog(group="Geometry"));
             Modelica.SIunits.SpecificEnthalpy initialEnthalpy[:]=
               {h_start, h_start}
              "Initial values of the Enthalpy for the spatialDistribuiton";

            Modelica.SIunits.Length x(start=0)
              "Spatial coordiante for spatialDistribution operator";
            Modelica.SIunits.SpecificEnthalpy h_start;
            input Modelica.SIunits.MassFlowRate m_flow;
            parameter Modelica.SIunits.Temp_K T_start;

          equation
            // start entalphy for the spatialDistribution, change this maybee to only do it once?
            h_start = Medium.specificEnthalpy_pTX(port_a.p,T_start,inStream(port_a.Xi_outflow));

            Volume = (diameter/2)*(diameter/2)*pi*length;
            delayTime = Volume*density/m_flow;
            density = Medium.density(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
            cp = Medium.specificHeatCapacityCp(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));

            der(x) = velocity;
            velocity = (length/delayTime);

            (h_b, h_a) =
                   spatialDistribution(
                     inStream(port_a.h_outflow),
                     inStream(port_b.h_outflow),
                     x/length,noEvent(velocity>=0),{0.0,1.0},initialEnthalpy);

          //    if velocity >= 0 then
          //     T_in = Medium.temperature(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));
          //    else
          //     T_in = Medium.temperature(Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow)));
          //    end if;

             if velocity >= 0 then
              T_in = Medium.temperature(Medium.setState_phX(port_b.p, h_a, inStream(port_a.Xi_outflow)));
             else
              T_in = Medium.temperature(Medium.setState_phX(port_a.p, h_b, inStream(port_b.Xi_outflow)));
            end if;

            (time_a, time_b) =
               spatialDistribution(time, time, x/length, noEvent(velocity>=0), {0.0,1.0},{0.0,0.0});

            if noEvent(velocity >= 0) then
              tau = time - time_b;
            else
              tau = time - time_a;
            end if;

            tauExp = tau/(R*density*pi*(diameter/2)^2*cp);
            exponential = exp(-tauExp);

            T_out = T_amb + (T_in - T_amb) * exponential;

            port_b.h_outflow = Medium.specificEnthalpy_pTX(port_b.p,T_out,inStream(port_a.Xi_outflow));
            port_a.h_outflow = Medium.specificEnthalpy_pTX(port_a.p,T_out,inStream(port_a.Xi_outflow));

            // No storage
            0 = port_a.m_flow + port_b.m_flow;

            // No pressure loss in this component
            port_b.p = port_a.p;

            port_b.Xi_outflow = inStream(port_a.Xi_outflow);
            port_b.C_outflow = inStream(port_a.C_outflow);

            port_a.Xi_outflow = inStream(port_b.Xi_outflow);
            port_a.C_outflow = inStream(port_b.C_outflow);
            annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                    Rectangle(
                    extent={{-100,60},{100,-60}},
                    fillColor={0,0,0},
                    fillPattern=FillPattern.Solid,
                    pattern=LinePattern.None,
                    lineColor={0,0,0}),
                  Text(
                    extent={{-145,-56},{155,-96}},
                    lineColor={0,140,72},
                    textString="%name"),
                  Rectangle(
                    extent={{-26,82},{28,-8}},
                    lineColor={255,0,0},
                    fillColor={255,0,0},
                    fillPattern=FillPattern.Solid),
                  Polygon(
                    points={{-40,80},{40,80},{0,126},{-40,80}},
                    lineColor={255,0,0},
                    fillColor={255,0,0},
                    fillPattern=FillPattern.Solid)}),                      Diagram(
                  coordinateSystem(preserveAspectRatio=false)),
                        Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
                  coordinateSystem(preserveAspectRatio=false)));
          end spatialDelayInput;
        end parts;
      end SpatialPipe;
    end Pipe;

    package Consumer

      model Substation
        "Consumer modl with PI control used for non BCVTB case"
      //   import DistictHeatingSimulation = DHNLibrary.DistrictHeating;
        //Works fine but there is a problem, it can´t handle small load, so when using ex. a ramp you have to give it a small offset like 0.01, something with the PID.
        Modelica.Fluid.Valves.ValveLinear valveLinear(
          redeclare package Medium = Medium,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=dp_nominal)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={34,0})));
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
              Medium)                                annotation (Placement(
              transformation(extent={{-112,-10},{-92,10}}), iconTransformation(extent=
                 {{-112,-10},{-92,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium =
              Medium)                                annotation (Placement(
              transformation(extent={{-110,-50},{-90,-30}}),
                                                          iconTransformation(extent={{-110,
                  -50},{-90,-30}})));
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
          annotation (__Dymola_choicesAllMatching=true);

        Modelica.Blocks.Interfaces.RealInput load
          annotation (Placement(transformation(extent={{142,60},{102,100}}),
              iconTransformation(extent={{122,48},{82,88}})));

        Modelica.Blocks.Continuous.LimPID PI(
          yMin=0,
          controllerType=Modelica.Blocks.Types.SimpleController.PI,
          yMax=1) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-24,60})));

        Modelica.Blocks.Interfaces.RealOutput HeatFlowHouse
          annotation (Placement(transformation(extent={{98,-16},{118,4}}),
              iconTransformation(extent={{98,-16},{118,4}})));
        parameter Modelica.SIunits.HeatFlowRate minimalPower = 0.03*abs(maximalPower)
          "Minimal heat flow rate that has to be used"
        annotation(Dialog(tab = "Advanced"));
        parameter Modelica.SIunits.HeatFlowRate maximalPower = 50000
          "The maximal heat flow rate that is able to be taken out from the house"
        annotation(Dialog(tab = "Advanced"));
        parameter Modelica.Media.Interfaces.PartialMedium.MassFlowRate m_flow_nominal
          "Nominal mass flowrate at full opening";
        parameter Modelica.SIunits.AbsolutePressure dp_nominal
          "Nominal pressure drop at full opening";

        Modelica.Blocks.Math.Product product
          annotation (Placement(transformation(extent={{-70,52},{-56,66}})));
        Modelica.Blocks.Sources.Constant kilo(k=1000)
          annotation (Placement(transformation(extent={{-112,22},{-92,42}})));
        Parts.minLoad
          minLoad(minLoad=minimalPower)
          annotation (Placement(transformation(extent={{-48,58},{-44,62}})));
        parameter Real k=1 "Gain" annotation (Dialog(tab="FirstOrder"));
        parameter Modelica.SIunits.Time T=1 "Time Constant"
          annotation (Dialog(tab="FirstOrder"));
        Modelica.Blocks.Continuous.FirstOrder firstOrder(k=k, T=T) annotation (
            Placement(transformation(
              extent={{-6,-6},{6,6}},
              rotation=90,
              origin={-20,20})));
        Parts.HeatExchangerFixedReturnT
          heatExchangerFixedReturnT(m_flow_nominal=m_flow_nominal,
          redeclare package Medium = Medium,
          T_return=T_return,
          pressure_out = port_b.p)
          annotation (Placement(transformation(extent={{-32,-12},{-16,4}})));
        parameter Modelica.SIunits.Temp_K T_return=303.15
          "Temperature at port_b for out-flowing fluid";
      equation

        connect(PI.y, valveLinear.opening) annotation (Line(points={{-13,60},{6,
                60},{6,8},{34,8}}, color={0,0,127}));
        connect(kilo.y, product.u2) annotation (Line(points={{-91,32},{-82,32},{-82,
                54.8},{-71.4,54.8}},
                               color={0,0,127}));
        connect(port_b, valveLinear.port_b) annotation (Line(points={{-100,-40},{-20,-40},
                {60,-40},{60,0},{44,0}}, color={0,127,255}));
        connect(product.y, minLoad.u) annotation (Line(points={{-55.3,59},{-51.65,59},
                {-51.65,60},{-47.6,60}}, color={0,0,127}));
        connect(minLoad.y, PI.u_s) annotation (Line(points={{-43.4,60},{-39.7,60},{-36,
                60}}, color={0,0,127}));
        connect(product.u1, load) annotation (Line(points={{-71.4,63.2},{-80.7,63.2},
                {-80.7,80},{122,80}}, color={0,0,127}));
        connect(HeatFlowHouse, PI.u_m) annotation (Line(points={{108,-6},{42,-6},{42,
                32},{-24,32},{-24,48}}, color={0,0,127}));
        connect(firstOrder.y, PI.u_m) annotation (Line(points={{-20,26.6},{-22,26.6},
                {-22,32},{-24,32},{-24,48}}, color={0,0,127}));
        connect(port_a, heatExchangerFixedReturnT.port_a)
          annotation (Line(points={{-102,0},{-32,0}}, color={0,127,255}));
        connect(heatExchangerFixedReturnT.port_b, valveLinear.port_a) annotation (
            Line(points={{-32,-8},{-40,-8},{-40,-18},{-4,-18},{-4,0},{24,0}}, color={0,
                127,255}));
        connect(heatExchangerFixedReturnT.HeatFlowHouse, firstOrder.u) annotation (
            Line(points={{-24,4.6},{-22,4.6},{-22,12.8},{-20,12.8}}, color={0,0,127}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
                extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}),
              Polygon(
                points={{0,200},{-160,90},{152,90},{0,200}},
                lineColor={0,0,0},
                fillColor={131,94,83},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-26,-24},{22,-100}},
                lineColor={118,91,72},
                fillColor={131,94,83},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-76,52},{-44,2}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                lineThickness=0.5),
              Rectangle(
                extent={{40,54},{72,4}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                lineThickness=0.5),
              Text(
                extent={{-143,-98},{157,-138}},
                lineColor={0,0,255},
                textString="%name")}));
      end Substation;

      model Consumer_BCVTB
        "Consumer model with a valve driven mass flow and return temperature inlet. Supply temperature is also output"
        import DistictHeatingSimulation = DHNLibrary.DistrictHeating;
        //Works fine but there is a problem, it can´t handle small load, so when using ex. a ramp you have to give it a small offset like 0.01, something with the PID.
        Modelica.Fluid.Valves.ValveLinear valveLinear(
          redeclare package Medium = Medium,
          m_flow_nominal=m_flow_nominal,
          dp_nominal=dp_nominal)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
              rotation=0,
              origin={34,0})));
        Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium
            = Medium)                                annotation (Placement(
              transformation(extent={{-112,-10},{-92,10}}), iconTransformation(extent=
                 {{-112,-10},{-92,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium
            = Medium)                                annotation (Placement(
              transformation(extent={{-110,-50},{-90,-30}}),
                                                          iconTransformation(extent={{-110,
                  -50},{-90,-30}})));
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
          annotation (__Dymola_choicesAllMatching=true);

        Modelica.Blocks.Interfaces.RealInput m_flow annotation (Placement(
              transformation(extent={{142,60},{102,100}}), iconTransformation(extent={{84,154},
                  {-4,242}})));

        parameter Modelica.SIunits.HeatFlowRate minimalPower = 0.03*abs(maximalPower)
          "Minimal heat flow rate that has to be used"
        annotation(Dialog(tab = "Advanced"));
        parameter Modelica.SIunits.HeatFlowRate maximalPower = 50000
          "The maximal heat flow rate that is able to be taken out from the house"
        annotation(Dialog(tab = "Advanced"));
        parameter Modelica.Media.Interfaces.PartialMedium.MassFlowRate m_flow_nominal
          "Nominal mass flowrate at full opening";
        parameter Modelica.SIunits.AbsolutePressure dp_nominal
          "Nominal pressure drop at full opening";

        parameter Real k=1 "Gain" annotation (Dialog(tab="FirstOrder"));
        parameter Modelica.SIunits.Time T=1 "Time Constant"
          annotation (Dialog(tab="FirstOrder"));
        parameter Modelica.SIunits.Temp_K T_return=303.15
          "Temperature at port_b for out-flowing fluid";
        Annex60.Fluid.HeatExchangers.HeaterCooler_T hea(
          redeclare package Medium = Medium,
          m_flow_nominal=5,
          dp_nominal=5000)
          annotation (Placement(transformation(extent={{0,-54},{-20,-34}})));
        Modelica.Blocks.Interfaces.RealInput T_re annotation (Placement(
              transformation(
              extent={{20,-20},{-20,20}},
              rotation=0,
              origin={124,-88}), iconTransformation(extent={{216,-158},{98,-40}})));
        Umrechnung umrechnung(m_flow_nominal=5)
          annotation (Placement(transformation(extent={{18,38},{38,58}})));
        Annex60.Fluid.Sensors.TemperatureTwoPort senTem(redeclare package
            Medium =
              Medium, m_flow_nominal=5)
          annotation (Placement(transformation(extent={{-46,-10},{-26,10}})));
        Modelica.Blocks.Interfaces.RealOutput T_Supply annotation (Placement(
              transformation(extent={{98,26},{118,46}}), iconTransformation(
                extent={{102,18},{172,88}})));
      equation

        connect(hea.port_b, port_b) annotation (Line(points={{-20,-44},{-60,-44},{-60,
                -40},{-100,-40}}, color={0,127,255}));
        connect(valveLinear.port_b, hea.port_a) annotation (Line(points={{44,0},{52,0},
                {52,2},{64,2},{64,-44},{0,-44}}, color={0,127,255}));
        connect(T_re, hea.TSet)
          annotation (Line(points={{124,-88},{2,-88},{2,-38}}, color={0,0,127}));
        connect(m_flow, umrechnung.m_flow) annotation (Line(points={{122,80},{74,80},
                {74,78},{28,78},{28,56}}, color={0,0,127}));
        connect(umrechnung.opening, valveLinear.opening) annotation (Line(points={{
                38.8,48},{46,48},{46,34},{34,34},{34,8}}, color={0,0,127}));
        connect(port_a, senTem.port_a) annotation (Line(points={{-102,0},{-74,0},
                {-46,0}}, color={0,127,255}));
        connect(valveLinear.port_a, senTem.port_b)
          annotation (Line(points={{24,0},{-26,0}}, color={0,127,255}));
        connect(senTem.T, T_Supply) annotation (Line(points={{-36,11},{32,11},{32,
                36},{108,36}}, color={0,0,127}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
                extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={255,85,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0}),
              Polygon(
                points={{0,200},{-160,90},{152,90},{0,200}},
                lineColor={0,0,0},
                fillColor={131,94,83},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-26,-24},{22,-100}},
                lineColor={118,91,72},
                fillColor={131,94,83},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-76,52},{-44,2}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                lineThickness=0.5),
              Rectangle(
                extent={{40,54},{72,4}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,0},
                lineThickness=0.5),
              Text(
                extent={{-143,-98},{157,-138}},
                lineColor={0,0,255},
                textString="%name")}));
      end Consumer_BCVTB;

      package Parts

        model HeatExchangerFixedReturnT

          Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium
              = Medium)
            annotation (Placement(transformation(extent={{-90,30},{-70,50}})));
          Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium
              = Medium)
            annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
          replaceable package Medium =
              Modelica.Media.Interfaces.PartialMedium
            annotation (__Dymola_choicesAllMatching=true);

          Modelica.SIunits.Temp_K T_in "Temperature at port_a for in-flowing fluid";
          parameter Modelica.SIunits.Temp_K T_return
            "Temperature at port_b for out-flowing fluid";

          parameter Modelica.SIunits.MassFlowRate m_flow_nominal
            "Nominal mass flow rate"
          annotation(Dialog(group = "Nominal condition"));
          Modelica.SIunits.MassFlowRate m_flow "Mass flow in the pipe";

          parameter Modelica.SIunits.MassFlowRate m_flow_small(min=0) = 1E-4*abs(m_flow_nominal)
            "Small mass flow rate for regularization of zero flow"
          annotation(Dialog(tab = "Advanced"));
          parameter Modelica.SIunits.MassFlowRate m_flow_min = 0.0002
            "Minimum mass flow rate for regularization of zero flow"
          annotation(Dialog(tab = "Advanced"));
          input Modelica.SIunits.AbsolutePressure pressure_out;

        public
          Modelica.Blocks.Interfaces.RealOutput HeatFlowHouse annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={0,86})));
        equation

          m_flow= noEvent(if port_a.m_flow > m_flow_min then port_a.m_flow else m_flow_small);

          T_in = Medium.temperature(Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)));

          port_a.h_outflow = Medium.specificEnthalpy_pTX(port_a.p,T_in,inStream(port_a.Xi_outflow));
          port_b.h_outflow = Medium.specificEnthalpy_pTX(pressure_out,T_return,inStream(port_a.Xi_outflow));

          HeatFlowHouse = (port_a.h_outflow - port_b.h_outflow)*port_a.m_flow;

          port_b.p = port_a.p;

          // No storage
          0 = port_a.m_flow + port_b.m_flow;

          port_b.Xi_outflow = inStream(port_a.Xi_outflow);
          port_b.C_outflow = inStream(port_a.C_outflow);

          port_a.Xi_outflow = inStream(port_b.Xi_outflow);
          port_a.C_outflow = inStream(port_b.C_outflow);

              annotation (
            Documentation(info="<html>
<p>A Simple Heat Exchanger used for customer in the district heating networt.</p>
<p>There is no delay between the incoming medium and outgoing, however there is a delay calculated to account for a higher</p><p>mass flow would increase the the temperature of the outgoing medium. </p>
<p>The HeatExchanger is only being tested without reverse flow (should not work with reverse)</p>
</html>"),  Diagram(coordinateSystem(extent={{-80,-80},{80,80}})),
            Icon(coordinateSystem(extent={{-80,-80},{80,80}}), graphics={
                                         Rectangle(
                  extent={{-80,80},{80,-80}},
                  pattern=LinePattern.None,
                  lineColor={0,0,0},
                  fillColor={255,255,0},
                  fillPattern=FillPattern.CrossDiag)}));
        end HeatExchangerFixedReturnT;

        model minLoad

          Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{16,-10},{36,10}})));
          Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
                  extent={{-42,-26},{-2,14}}), iconTransformation(extent={{-30,-14},{-2,
                    14}})));
         parameter Real minLoad = 5;
        equation
          y= noEvent(if u > minLoad then u else minLoad);
          annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-20,-20},
                    {20,20}}), graphics={Rectangle(
                  extent={{-22,22},{22,-20}},
                  lineColor={28,108,200},
                  fillColor={170,255,255},
                  fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                  preserveAspectRatio=false, extent={{-20,-20},{20,20}})));
        end minLoad;
      end Parts;
    end Consumer;
  end Components;

  package Experiments

    model BCVTBblock
      Buildings.Utilities.IO.BCVTB.BCVTB bcvtb(
        xmlFileName="socket.cfg",
        timeStep=60,
        final nDblRea=34,
        final nDblWri=19,
        uStart={60,30,30,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59})
        annotation (Placement(transformation(extent={{-70,-2},{-50,18}})));
      Modelica.Blocks.Routing.Multiplex2 write(n1=3, n2=16)
        annotation (Placement(transformation(extent={{210,0},{230,20}})));
      Modelica.Blocks.Routing.DeMultiplex2 read(n1=2, n2=32)
        annotation (Placement(transformation(extent={{-8,-2},{12,18}})));
    equation
      connect(bcvtb.yR,read. u) annotation (Line(
          points={{-49,8},{-10,8}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(write.y,bcvtb. uR) annotation (Line(
          points={{231,10},{240,10},{240,-82},{-80,-82},{-80,8},{-72,8}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{260,100}})), Diagram(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{260,100}})));
    end BCVTBblock;

    model GrazNetworkModel1
      "24 hour network simulation driven my mass flow source with step increase in mass flow and supply temperature "

      Components.Consumer.Substation K(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{2,76},{14,90}})));
      Components.Consumer.Substation J(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-46,76},{-60,90}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L1_2(
        redeclare package Medium = Medium,
        diameter=0.2065,
        InsThickness=0.04975,
        m_flow_nominal=100,
        length=2000,
        T_start=333.15)
        annotation (Placement(transformation(extent={{6,-12},{26,8}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L3(
        length=102,
        redeclare package Medium = Medium,
        diameter=0.2065,
        InsThickness=0.04975,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-92,-14},{-72,6}})));
      Components.Consumer.Substation L(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-116,-46},{-130,-32}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L4(
        length=107,
        redeclare package Medium = Medium,
        diameter=0.2065,
        InsThickness=0.04975,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-194,-14},{-174,6}})));
      Components.Consumer.Substation B(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-230,74},{-244,88}})));
      Components.Consumer.Substation C(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-214,-54},{-228,-40}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L5(
        length=82,
        redeclare package Medium = Medium,
        diameter=0.1593,
        InsThickness=0.04135,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-268,-14},{-248,6}})));
      Components.Consumer.Substation A(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-252,52},{-238,66}})));
      Components.Consumer.Substation E(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-308,104},{-322,118}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L13(
        length=302,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-9,-8},{9,8}},
            rotation=-90,
            origin={-20,39})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L6(
        length=68,
        redeclare package Medium = Medium,
        diameter=0.1593,
        InsThickness=0.04135,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-340,-16},{-320,4}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L15(
        length=276,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-292,33})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L17(
        length=85,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-8,7},{8,-7}},
            rotation=-90,
            origin={-291,80})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L18(
        length=185,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-8,7},{8,-7}},
            rotation=-90,
            origin={-293,146})));
      Components.Consumer.Substation D(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-276,170},{-262,184}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L7(
        length=77,
        redeclare package Medium = Medium,
        diameter=0.1317,
        InsThickness=0.04305,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-408,-16},{-390,4}})));
      Components.Consumer.Substation I(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-434,-80},{-448,-66}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L19a_b(
        length=356,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,7},{7,-7}},
            rotation=-90,
            origin={-419,-49})));
      Components.Consumer.Substation H(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-360,-52},{-378,-34}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L8(
        length=66,
        redeclare package Medium = Medium,
        diameter=0.1317,
        InsThickness=0.04305,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-460,-16},{-442,4}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L9(
        length=55,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-522,-16},{-504,4}})));
      Components.Consumer.Substation F(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-480,24},{-494,38}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L10(
        length=79,
        redeclare package Medium = Medium,
        diameter=0.0825,
        InsThickness=0.03555,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-558,-18},{-540,2}})));
      Components.Consumer.Substation G(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-536,-50},{-550,-36}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L11(
        length=72,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-596,-18},{-578,2}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L20(
        length=81,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-566,27})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L21(
        length=184,
        redeclare package Medium = Medium,
        diameter=0.0372,
        InsThickness=0.0338,
        m_flow_nominal=100,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-566,69})));
      Components.Consumer.Substation N(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-590,42},{-604,56}})));
      Components.Consumer.Substation M(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-582,94},{-596,108}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L12(
        length=83,
        redeclare package Medium = Medium,
        diameter=0.0545,
        InsThickness=0.03235,
        m_flow_nominal=100,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-646,-18},{-628,2}})));
      Components.Consumer.Substation P(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-612,-56},{-626,-42}})));
      Components.Consumer.Substation O(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-668,10},{-654,24}})));
      parameter Modelica.SIunits.ThermalConductivity lambdaG=2.5
        "Thermal conductivity of ground layers";
      Annex60.Fluid.Sources.Boundary_pT Sink(
        nPorts=1,
        redeclare package Medium = Medium,
        p=100000,
        T=313.15) annotation (Placement(transformation(extent={{88,10},{68,30}})));
      replaceable package Medium =  Annex60.Media.Water;
      Modelica.Blocks.Sources.Constant
                                   Load(k=1000)
        "Constant load to each consumer"
        annotation (Placement(transformation(extent={{50,130},{30,150}})));
      Modelica.Blocks.Sources.Step TempStepChange(
        height=10,
        offset=273 + 60,
        startTime=8000)
        annotation (Placement(transformation(extent={{140,58},{160,78}})));
      Annex60.Fluid.Sources.MassFlowSource_T MassFlowSource(
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1,
        redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{102,-22},{82,-2}})));
      Modelica.Blocks.Sources.Step MassFlowStepUp(
        height=5,
        offset=100,
        startTime=8000)
        annotation (Placement(transformation(extent={{84,70},{104,90}})));
    equation
      connect(L5.port_a2, L4.port_b2) annotation (Line(points={{-248,-10},{-222,
              -10},{-194,-10}}, color={0,127,255}));
      connect(L5.port_b1, L4.port_a1) annotation (Line(points={{-248,0},{-218,0},
              {-194,0}}, color={0,127,255}));
      connect(L3.port_a2, L1_2.port_b2) annotation (Line(points={{-72,-10},{10,
              -10},{6,-8}},   color={0,127,255}));
      connect(L3.port_b1, L1_2.port_a1)
        annotation (Line(points={{-72,0},{10,0},{6,2}},   color={0,127,255}));
      connect(L13.port_b1, L1_2.port_a1) annotation (Line(points={{-16.8,30},{
              -16.8,2},{6,2}},  color={0,127,255}));
      connect(L1_2.port_b2, L13.port_a2) annotation (Line(points={{6,-8},{-24,
              -8},{-24.8,-8},{-24.8,30}},   color={0,127,255}));
      connect(L6.port_b1, L5.port_a1) annotation (Line(points={{-320,-2},{-316,
              0},{-268,0}},
                         color={0,127,255}));
      connect(L6.port_a2, L5.port_b2) annotation (Line(points={{-320,-12},{-316,
              -10},{-268,-10}}, color={0,127,255}));
      connect(L15.port_a2, L5.port_b2) annotation (Line(points={{-287.2,26},{
              -284,26},{-284,-10},{-268,-10}},
                                          color={0,127,255}));
      connect(L15.port_b1, L5.port_a1) annotation (Line(points={{-295.2,26},{
              -294,26},{-294,0},{-268,0}},
                                      color={0,127,255}));
      connect(L17.port_b1, L15.port_a1) annotation (Line(points={{-293.8,72},{
              -294,72},{-294,64},{-295.2,64},{-295.2,40}}, color={0,127,255}));
      connect(L17.port_a2, L15.port_b2) annotation (Line(points={{-286.8,72},{
              -286.8,56},{-287.2,56},{-287.2,40}}, color={0,127,255}));
      connect(L17.port_b2, L18.port_a2) annotation (Line(points={{-286.8,88},{
              -288,88},{-288,138},{-288.8,138}}, color={0,127,255}));
      connect(L17.port_a1, L18.port_b1) annotation (Line(points={{-293.8,88},{
              -296,88},{-296,138},{-295.8,138}}, color={0,127,255}));
      connect(L3.port_a1, L4.port_b1)
        annotation (Line(points={{-92,0},{-86,0},{-174,0}}, color={0,127,255}));
      connect(L4.port_a2, L3.port_b2)
        annotation (Line(points={{-174,-10},{-92,-10}}, color={0,127,255}));
      connect(L7.port_b1, L6.port_a1)
        annotation (Line(points={{-390,-2},{-374,-2},{-374,-4},{-362,-4},{-362,
              -2},{-340,-2}},                        color={0,127,255}));
      connect(L7.port_a2, L6.port_b2)
        annotation (Line(points={{-390,-12},{-340,-12}}, color={0,127,255}));
      connect(L8.port_b1, L7.port_a1)
        annotation (Line(points={{-442,-2},{-430,-2},{-430,-4},{-422,-4},{-422,
              -2},{-408,-2}},                        color={0,127,255}));
      connect(L8.port_a2, L7.port_b2) annotation (Line(points={{-442,-12},{-422,
              -12},{-408,-12}}, color={0,127,255}));
      connect(L19a_b.port_a1, L7.port_b2) annotation (Line(points={{-421.8,-42},
              {-421.8,-12},{-408,-12}},color={0,127,255}));
      connect(L19a_b.port_b2, L7.port_a1) annotation (Line(points={{-414.8,-42},
              {-414,-42},{-414,-2},{-408,-2}},
                                            color={0,127,255}));
      connect(I.port_a, L19a_b.port_b1) annotation (Line(points={{-433.86,-73},
              {-421.8,-73},{-421.8,-56}},color={0,127,255}));
      connect(I.port_b, L19a_b.port_a2) annotation (Line(points={{-434,-75.8},{
              -424,-75.8},{-424,-74},{-414.8,-74},{-414.8,-56}}, color={0,127,255}));
      connect(L9.port_b1, L8.port_a1)
        annotation (Line(points={{-504,-2},{-490,-2},{-490,-4},{-478,-4},{-478,
              -2},{-460,-2}},                        color={0,127,255}));
      connect(L9.port_a2, L8.port_b2)
        annotation (Line(points={{-504,-12},{-460,-12}}, color={0,127,255}));
      connect(L10.port_b1, L9.port_a1) annotation (Line(points={{-540,-4},{-536,
              -2},{-522,-2}},
                         color={0,127,255}));
      connect(L10.port_a2, L9.port_b2) annotation (Line(points={{-540,-14},{
              -536,-12},{-522,-12}},
                                color={0,127,255}));
      connect(L11.port_b1, L10.port_a1) annotation (Line(points={{-578,-4},{
              -565,-4},{-558,-4}},
                              color={0,127,255}));
      connect(L11.port_a2, L10.port_b2) annotation (Line(points={{-578,-14},{
              -565,-14},{-558,-14}},
                                color={0,127,255}));
      connect(L20.port_a2, L10.port_b2) annotation (Line(points={{-561.2,20},{
              -560,20},{-560,-14},{-558,-14}}, color={0,127,255}));
      connect(L20.port_b1, L10.port_a1) annotation (Line(points={{-569.2,20},{
              -570,20},{-570,-4},{-558,-4}}, color={0,127,255}));
      connect(L21.port_a2, L20.port_b2) annotation (Line(points={{-561.2,62},{
              -562,62},{-562,34},{-561.2,34}}, color={0,127,255}));
      connect(L20.port_a1, L21.port_b1) annotation (Line(points={{-569.2,34},{
              -570.8,34},{-570.8,62},{-569.2,62}},      color={0,127,255}));
      connect(N.port_a, L20.port_b2) annotation (Line(points={{-589.86,49},{-562,49},
              {-562,34},{-561.2,34}},     color={0,127,255}));
      connect(N.port_b, L21.port_b1) annotation (Line(points={{-590,46.2},{-582,
              46.2},{-582,46},{-570,46},{-569.2,46},{-569.2,62}}, color={0,127,
              255}));
      connect(M.port_b, L21.port_a1) annotation (Line(points={{-582,98.2},{-570,
              98.2},{-570,92},{-570,76},{-569.2,76}}, color={0,127,255}));
      connect(M.port_a, L21.port_b2) annotation (Line(points={{-581.86,101},{-561.2,
              101},{-561.2,76}},        color={0,127,255}));
      connect(C.port_a, L4.port_b2) annotation (Line(points={{-213.86,-47},{
              -206,-45},{-206,-10},{-194,-10}},
                                           color={0,127,255}));
      connect(C.port_b, L4.port_a1) annotation (Line(points={{-214,-49.8},{-208,
              -49.8},{-208,-48},{-202,-48},{-202,0},{-194,0}}, color={0,127,255}));
      connect(L12.port_b1, L11.port_a1) annotation (Line(points={{-628,-4},{
              -610,-4},{-596,-4}},
                              color={0,127,255}));
      connect(L12.port_a2, L11.port_b2) annotation (Line(points={{-628,-14},{
              -610,-14},{-596,-14}},
                                color={0,127,255}));
      connect(L.port_a, L3.port_b2) annotation (Line(points={{-115.86,-39},{
              -110,-37},{-110,-10},{-92,-10}},
                                          color={0,127,255}));
      connect(L.port_b, L4.port_b1) annotation (Line(points={{-116,-41.8},{-110,
              -41.8},{-110,-40},{-106,-40},{-106,0},{-174,0}}, color={0,127,255}));
      connect(H.port_a, L6.port_b2) annotation (Line(points={{-359.82,-43},{
              -352,-43},{-352,-12},{-340,-12}},
                                           color={0,127,255}));
      connect(H.port_b, L6.port_a1) annotation (Line(points={{-360,-46.6},{-352,
              -46.6},{-352,-48},{-346,-48},{-346,-2},{-340,-2}},
                                                               color={0,127,255}));
      connect(G.port_b, L9.port_a1) annotation (Line(points={{-536,-45.8},{-532,
              -45.8},{-532,-46},{-524,-46},{-524,-2},{-522,-2}},
                                                               color={0,127,255}));
      connect(G.port_a, L9.port_b2) annotation (Line(points={{-535.86,-43},{
              -528,-43},{-528,-12},{-522,-12}},
                                           color={0,127,255}));
      connect(P.port_b, L11.port_a1) annotation (Line(points={{-612,-51.8},{
              -606,-51.8},{-606,-52},{-598,-52},{-598,-4},{-596,-4}},
                                                                 color={0,127,255}));
      connect(P.port_a, L11.port_b2) annotation (Line(points={{-611.86,-49},{
              -604,-49},{-604,-14},{-596,-14}},
                                           color={0,127,255}));
      connect(L1_2.port_b1, Sink.ports[1]) annotation (Line(points={{26,2},{44,
              2},{44,20},{68,20}},
                                color={0,127,255}));
      connect(F.port_b, L8.port_a1) annotation (Line(points={{-480,28.2},{-478,
              28.2},{-478,28},{-474,28},{-474,-2},{-460,-2}},
                                                            color={0,127,255}));
      connect(F.port_a, L8.port_b2) annotation (Line(points={{-479.86,31},{-466,
              31},{-466,-12},{-460,-12}}, color={0,127,255}));
      connect(B.port_b, L4.port_a1) annotation (Line(points={{-230,78.2},{-226,
              78.2},{-226,74},{-222,74},{-222,0},{-194,0}}, color={0,127,255}));
      connect(B.port_a, L4.port_b2) annotation (Line(points={{-229.86,81},{-214,
              81},{-214,-10},{-194,-10}}, color={0,127,255}));
      connect(E.port_b, L18.port_b1) annotation (Line(points={{-308,108.2},{
              -302,108.2},{-302,106},{-296,106},{-296,138},{-295.8,138}},
                                                                     color={0,127,
              255}));
      connect(E.port_a, L18.port_a2) annotation (Line(points={{-307.86,111},{
              -288,111},{-288,138},{-288.8,138}},
                                             color={0,127,255}));
      connect(D.port_b, L18.port_a1) annotation (Line(points={{-276,174.2},{
              -286,174.2},{-286,174},{-295.8,174},{-295.8,154}},
                                                            color={0,127,255}));
      connect(L18.port_b2, D.port_a) annotation (Line(points={{-288.8,154},{
              -288,154},{-288,180},{-288,177},{-276.14,177}},
                                                         color={0,127,255}));
      connect(L13.port_a1, J.port_b) annotation (Line(points={{-16.8,48},{-16,
              48},{-16,82},{-26.8,82},{-26.8,80.2},{-46,80.2}},
                                                            color={0,127,255}));
      connect(L13.port_b2, K.port_a) annotation (Line(points={{-24.8,48},{-28,
              48},{-28,83},{-18,83},{1.88,83}},
                                            color={0,127,255}));
      connect(J.port_a, K.port_a) annotation (Line(points={{-45.86,83},{-28,85},{1.88,
              83}},      color={0,127,255}));
      connect(K.port_b, J.port_b) annotation (Line(points={{2,80.2},{-8,82.2},{-8,78},
              {-16,78},{-16,82},{-26.8,82},{-26.8,82.2},{-46,80.2}},     color={0,
              127,255}));
      connect(O.port_b, L12.port_a1) annotation (Line(points={{-668,14.2},{-672,
              14.2},{-672,14},{-676,14},{-676,-4},{-646,-4}}, color={0,127,255}));
      connect(O.port_a, L12.port_b2) annotation (Line(points={{-668.14,17},{
              -680,17},{-680,-14},{-646,-14}},
                                          color={0,127,255}));
      connect(Load.y, K.load) annotation (Line(points={{29,140},{22,140},{22,
              87.76},{14.12,87.76}}, color={0,0,127}));
      connect(Load.y, J.load) annotation (Line(points={{29,140},{-14,140},{-14,
              87.76},{-60.14,87.76}}, color={0,0,127}));
      connect(Load.y, L.load) annotation (Line(points={{29,140},{-50,140},{-50,
              -34.24},{-130.14,-34.24}}, color={0,0,127}));
      connect(Load.y, C.load) annotation (Line(points={{29,140},{-100,140},{
              -100,-42.24},{-228.14,-42.24}}, color={0,0,127}));
      connect(Load.y, B.load) annotation (Line(points={{29,140},{-106,140},{
              -106,85.76},{-244.14,85.76}}, color={0,0,127}));
      connect(Load.y, A.load) annotation (Line(points={{29,140},{-114,140},{
              -114,63.76},{-237.86,63.76}}, color={0,0,127}));
      connect(Load.y, D.load) annotation (Line(points={{29,140},{-116,140},{
              -116,181.76},{-261.86,181.76}}, color={0,0,127}));
      connect(Load.y, E.load) annotation (Line(points={{29,140},{-146,140},{
              -146,115.76},{-322.14,115.76}}, color={0,0,127}));
      connect(Load.y, H.load) annotation (Line(points={{29,140},{-174,140},{
              -174,-36.88},{-378.18,-36.88}}, color={0,0,127}));
      connect(Load.y, F.load) annotation (Line(points={{29,140},{-232,140},{
              -232,35.76},{-494.14,35.76}}, color={0,0,127}));
      connect(Load.y, I.load) annotation (Line(points={{29,140},{-210,140},{
              -210,-68.24},{-448.14,-68.24}}, color={0,0,127}));
      connect(Load.y, M.load) annotation (Line(points={{29,140},{-282.5,140},{
              -282.5,105.76},{-596.14,105.76}}, color={0,0,127}));
      connect(Load.y, N.load) annotation (Line(points={{29,140},{-286,140},{
              -286,53.76},{-604.14,53.76}}, color={0,0,127}));
      connect(Load.y, O.load) annotation (Line(points={{29,140},{-312,140},{
              -312,21.76},{-653.86,21.76}}, color={0,0,127}));
      connect(Load.y, P.load) annotation (Line(points={{29,140},{-302,140},{
              -302,-44.24},{-626.14,-44.24}}, color={0,0,127}));
      connect(Load.y, G.load) annotation (Line(points={{29,140},{-260,140},{
              -260,-38.24},{-550.14,-38.24}}, color={0,0,127}));
      connect(A.port_a, L15.port_b2) annotation (Line(points={{-252.14,59},{
              -284,59},{-284,56},{-287.2,56},{-287.2,40}}, color={0,127,255}));
      connect(A.port_b, L15.port_a1) annotation (Line(points={{-252,56.2},{-270,
              56.2},{-270,50},{-295.2,50},{-295.2,40}}, color={0,127,255}));
      connect(MassFlowSource.ports[1], L1_2.port_a2) annotation (Line(points={{
              82,-12},{58,-12},{58,-8},{26,-8}}, color={0,127,255}));
      connect(TempStepChange.y, MassFlowSource.T_in) annotation (Line(points={{
              161,68},{136,68},{136,-8},{104,-8}}, color={0,0,127}));
      connect(MassFlowStepUp.y, MassFlowSource.m_flow_in) annotation (Line(
            points={{105,80},{118,80},{118,82},{128,82},{128,-4},{102,-4}},
            color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-720,
                -180},{180,240}})), Diagram(coordinateSystem(preserveAspectRatio=
                false, extent={{-720,-180},{180,240}})));
    end GrazNetworkModel1;

    model GrazNetworkModel2 "Graz network connected to BCVTB for cosimulation"

      Components.Consumer.Consumer_BCVTB K(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{2,76},{14,90}})));
      Components.Consumer.Consumer_BCVTB J(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-46,76},{-60,90}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L1_2(
        redeclare package Medium = Medium,
        m_flow_nominal=50,
        diameter=0.2065,
        InsThickness=0.04975,
        length=192,
        T_start=333.15)
        annotation (Placement(transformation(extent={{12,-14},{32,6}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L3(
        length=102,
        redeclare package Medium = Medium,
        diameter=0.2065,
        InsThickness=0.04975,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-88,-14},{-68,6}})));
      Components.Consumer.Consumer_BCVTB L(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-116,-46},{-130,-32}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L4(
        m_flow_nominal=1,
        length=107,
        redeclare package Medium = Medium,
        diameter=0.2065,
        InsThickness=0.04975,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-190,-14},{-170,6}})));
      Components.Consumer.Consumer_BCVTB B(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-172,124},{-186,138}})));
      Components.Consumer.Consumer_BCVTB C(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-214,-54},{-228,-40}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L5(
        length=82,
        redeclare package Medium = Medium,
        diameter=0.1593,
        InsThickness=0.04135,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-262,-14},{-242,6}})));
      Components.Consumer.Consumer_BCVTB A(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-238,52},{-252,66}})));
      Components.Consumer.Consumer_BCVTB E(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-308,104},{-322,118}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L13(
        length=302,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=30,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-9,-8},{9,8}},
            rotation=-90,
            origin={-22,39})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L6(
        length=68,
        redeclare package Medium = Medium,
        diameter=0.1593,
        InsThickness=0.04135,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-336,-16},{-316,4}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L15(
        length=276,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-288,33})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L17(
        length=85,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-8,7},{8,-7}},
            rotation=-90,
            origin={-291,80})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L18(
        length=185,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-8,7},{8,-7}},
            rotation=-90,
            origin={-293,138})));
      Components.Consumer.Consumer_BCVTB D(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-276,170},{-262,184}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L7(
        length=77,
        redeclare package Medium = Medium,
        diameter=0.1317,
        InsThickness=0.04305,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-404,-16},{-386,4}})));
      Components.Consumer.Consumer_BCVTB I(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-434,-80},{-448,-66}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L19a_b(
        length=356,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,7},{7,-7}},
            rotation=-90,
            origin={-421,-47})));
      Components.Consumer.Consumer_BCVTB H(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-360,-52},{-378,-34}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L8(
        length=66,
        redeclare package Medium = Medium,
        diameter=0.1317,
        InsThickness=0.04305,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-456,-16},{-438,4}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L9(
        length=55,
        redeclare package Medium = Medium,
        diameter=0.1071,
        InsThickness=0.04285,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-518,-16},{-500,4}})));
      Components.Consumer.Consumer_BCVTB F(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-480,24},{-494,38}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L10(
        length=79,
        redeclare package Medium = Medium,
        diameter=0.0825,
        InsThickness=0.03555,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-554,-18},{-536,2}})));
      Components.Consumer.Consumer_BCVTB G(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-536,-50},{-550,-36}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L11(
        length=72,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-594,-18},{-576,2}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L20(
        length=81,
        redeclare package Medium = Medium,
        diameter=0.0703,
        InsThickness=0.03195,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-566,27})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L21(
        length=184,
        redeclare package Medium = Medium,
        diameter=0.0372,
        InsThickness=0.0338,
        m_flow_nominal=20,
        T_start=333.15) annotation (Placement(transformation(
            extent={{-7,8},{7,-8}},
            rotation=-90,
            origin={-566,69})));
      Components.Consumer.Consumer_BCVTB N(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-590,42},{-604,56}})));
      Components.Consumer.Consumer_BCVTB M(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-582,94},{-596,108}})));
      Components.Pipe.SpatialPipe.SpatialPipeVolume_DHN L12(
        length=83,
        redeclare package Medium = Medium,
        diameter=0.0545,
        InsThickness=0.03235,
        m_flow_nominal=50,
        T_start=333.15)
        annotation (Placement(transformation(extent={{-644,-18},{-626,2}})));
      Components.Consumer.Consumer_BCVTB P(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-612,-56},{-626,-42}})));
      Components.Consumer.Consumer_BCVTB O(
        m_flow_nominal=1,
        redeclare package Medium = Medium,
        dp_nominal=300000,
        T_return=316.15)
        annotation (Placement(transformation(extent={{-668,10},{-654,24}})));
      parameter Modelica.SIunits.ThermalConductivity lambdaG=2.5
        "Thermal conductivity of ground layers";
      Annex60.Fluid.Sources.Boundary_pT Sink(
        redeclare package Medium = Medium,
        p=100000,
        T=313.15,
        nPorts=1) annotation (Placement(transformation(extent={{220,2},{200,22}})));
      replaceable package Medium =  Annex60.Media.Water;
      Annex60.Fluid.Sources.MassFlowSource_T MassFlowSource(
        use_m_flow_in=true,
        use_T_in=true,
        nPorts=1,
        redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{102,-22},{82,-2}})));
      Buildings.Utilities.IO.BCVTB.BCVTB BCTVB(
        xmlFileName="socket.cfg",
        timeStep=60,
        final nDblRea=34,
        final nDblWri=18,
        uStart={60,30,30,59,59,59,59,59,59,59,59,59,59,59,59,59,59,59})
        annotation (Placement(transformation(extent={{-1080,-342},{-864,-162}})));
      Modelica.Blocks.Routing.Multiplex2 write(      n2=16, n1=2)
        annotation (Placement(transformation(extent={{-60,164},{68,292}})));
      Modelica.Blocks.Routing.DeMultiplex2 read(n1=2, n2=32)
        annotation (Placement(transformation(extent={{-738,178},{-634,282}})));
      Annex60.Fluid.Sensors.MassFlowRate MassFlowSens(redeclare package Medium =
            Medium)
        annotation (Placement(transformation(extent={{42,6},{62,26}})));
      Annex60.Fluid.Sensors.TemperatureTwoPort TempSens(redeclare package
          Medium = Medium, m_flow_nominal=50)
        annotation (Placement(transformation(extent={{142,8},{162,28}})));
    equation
      connect(L5.port_a2, L4.port_b2) annotation (Line(points={{-242,-10},{-222,
              -10},{-190,-10}}, color={0,127,255}));
      connect(L5.port_b1, L4.port_a1) annotation (Line(points={{-242,0},{-218,0},
              {-190,0}}, color={0,127,255}));
      connect(L3.port_a2, L1_2.port_b2) annotation (Line(points={{-68,-10},{-66,
              -10},{12,-10}}, color={0,127,255}));
      connect(L3.port_b1, L1_2.port_a1)
        annotation (Line(points={{-68,0},{-66,0},{12,0}}, color={0,127,255}));
      connect(L13.port_b1, L1_2.port_a1) annotation (Line(points={{-18.8,30},{
              -18.8,0},{12,0}}, color={0,127,255}));
      connect(L1_2.port_b2, L13.port_a2) annotation (Line(points={{12,-10},{-24,
              -10},{-26.8,-10},{-26.8,30}}, color={0,127,255}));
      connect(L6.port_b1, L5.port_a1) annotation (Line(points={{-316,-2},{-316,
              0},{-262,0}},
                         color={0,127,255}));
      connect(L6.port_a2, L5.port_b2) annotation (Line(points={{-316,-12},{-316,
              -10},{-262,-10}}, color={0,127,255}));
      connect(L15.port_a2, L5.port_b2) annotation (Line(points={{-283.2,26},{
              -284,26},{-284,-10},{-262,-10}},
                                          color={0,127,255}));
      connect(L15.port_b1, L5.port_a1) annotation (Line(points={{-291.2,26},{
              -294,26},{-294,0},{-262,0}},
                                      color={0,127,255}));
      connect(L17.port_b1, L15.port_a1) annotation (Line(points={{-293.8,72},{
              -294,72},{-294,64},{-291.2,64},{-291.2,40}}, color={0,127,255}));
      connect(L17.port_a2, L15.port_b2) annotation (Line(points={{-286.8,72},{
              -286.8,56},{-283.2,56},{-283.2,40}}, color={0,127,255}));
      connect(L17.port_b2, L18.port_a2) annotation (Line(points={{-286.8,88},{
              -288,88},{-288,130},{-288.8,130}}, color={0,127,255}));
      connect(L17.port_a1, L18.port_b1) annotation (Line(points={{-293.8,88},{
              -296,88},{-296,130},{-295.8,130}}, color={0,127,255}));
      connect(L3.port_a1, L4.port_b1)
        annotation (Line(points={{-88,0},{-86,0},{-170,0}}, color={0,127,255}));
      connect(L4.port_a2, L3.port_b2)
        annotation (Line(points={{-170,-10},{-88,-10}}, color={0,127,255}));
      connect(L7.port_b1, L6.port_a1)
        annotation (Line(points={{-386,-2},{-374,-2},{-374,-4},{-362,-4},{-362,
              -2},{-336,-2}},                        color={0,127,255}));
      connect(L7.port_a2, L6.port_b2)
        annotation (Line(points={{-386,-12},{-336,-12}}, color={0,127,255}));
      connect(L8.port_b1, L7.port_a1)
        annotation (Line(points={{-438,-2},{-430,-2},{-430,-4},{-422,-4},{-422,
              -2},{-404,-2}},                        color={0,127,255}));
      connect(L8.port_a2, L7.port_b2) annotation (Line(points={{-438,-12},{-422,
              -12},{-404,-12}}, color={0,127,255}));
      connect(L19a_b.port_a1, L7.port_b2) annotation (Line(points={{-423.8,-40},
              {-423.8,-12},{-404,-12}},color={0,127,255}));
      connect(L19a_b.port_b2, L7.port_a1) annotation (Line(points={{-416.8,-40},
              {-414,-40},{-414,-2},{-404,-2}},
                                            color={0,127,255}));
      connect(I.port_a, L19a_b.port_b1) annotation (Line(points={{-433.86,-73},
              {-423.8,-73},{-423.8,-54}},color={0,127,255}));
      connect(I.port_b, L19a_b.port_a2) annotation (Line(points={{-434,-75.8},{-424,
              -75.8},{-424,-74},{-416.8,-74},{-416.8,-54}},      color={0,127,255}));
      connect(L9.port_b1, L8.port_a1)
        annotation (Line(points={{-500,-2},{-490,-2},{-490,-4},{-478,-4},{-478,
              -2},{-456,-2}},                        color={0,127,255}));
      connect(L9.port_a2, L8.port_b2)
        annotation (Line(points={{-500,-12},{-456,-12}}, color={0,127,255}));
      connect(L10.port_b1, L9.port_a1) annotation (Line(points={{-536,-4},{-536,
              -2},{-518,-2}},
                         color={0,127,255}));
      connect(L10.port_a2, L9.port_b2) annotation (Line(points={{-536,-14},{-536,
              -12},{-518,-12}}, color={0,127,255}));
      connect(L11.port_b1, L10.port_a1) annotation (Line(points={{-576,-4},{
              -565,-4},{-554,-4}},
                              color={0,127,255}));
      connect(L11.port_a2, L10.port_b2) annotation (Line(points={{-576,-14},{-565,
              -14},{-554,-14}}, color={0,127,255}));
      connect(L20.port_a2, L10.port_b2) annotation (Line(points={{-561.2,20},{
              -560,20},{-560,-14},{-554,-14}}, color={0,127,255}));
      connect(L20.port_b1, L10.port_a1) annotation (Line(points={{-569.2,20},{
              -570,20},{-570,-4},{-554,-4}}, color={0,127,255}));
      connect(L21.port_a2, L20.port_b2) annotation (Line(points={{-561.2,62},{
              -562,62},{-562,34},{-561.2,34}}, color={0,127,255}));
      connect(L20.port_a1, L21.port_b1) annotation (Line(points={{-569.2,34},{
              -570.8,34},{-570.8,62},{-569.2,62}},      color={0,127,255}));
      connect(N.port_a, L20.port_b2) annotation (Line(points={{-589.86,49},{-562,49},
              {-562,34},{-561.2,34}},     color={0,127,255}));
      connect(N.port_b, L21.port_b1) annotation (Line(points={{-590,46.2},{-582,
              46.2},{-582,46},{-570,46},{-569.2,46},{-569.2,62}}, color={0,127,
              255}));
      connect(M.port_b, L21.port_a1) annotation (Line(points={{-582,98.2},{-570,
              98.2},{-570,92},{-570,76},{-569.2,76}}, color={0,127,255}));
      connect(M.port_a, L21.port_b2) annotation (Line(points={{-581.86,101},{-561.2,
              101},{-561.2,76}},        color={0,127,255}));
      connect(C.port_a, L4.port_b2) annotation (Line(points={{-213.86,-47},{
              -206,-45},{-206,-10},{-190,-10}},
                                           color={0,127,255}));
      connect(C.port_b, L4.port_a1) annotation (Line(points={{-214,-49.8},{-208,
              -49.8},{-208,-48},{-202,-48},{-202,0},{-190,0}}, color={0,127,255}));
      connect(L12.port_b1, L11.port_a1) annotation (Line(points={{-626,-4},{
              -610,-4},{-594,-4}},
                              color={0,127,255}));
      connect(L12.port_a2, L11.port_b2) annotation (Line(points={{-626,-14},{-610,
              -14},{-594,-14}}, color={0,127,255}));
      connect(L.port_a, L3.port_b2) annotation (Line(points={{-115.86,-39},{
              -110,-37},{-110,-10},{-88,-10}},
                                          color={0,127,255}));
      connect(L.port_b, L4.port_b1) annotation (Line(points={{-116,-41.8},{-110,
              -41.8},{-110,-40},{-106,-40},{-106,0},{-170,0}}, color={0,127,255}));
      connect(H.port_a, L6.port_b2) annotation (Line(points={{-359.82,-43},{-352,-43},
              {-352,-12},{-336,-12}},      color={0,127,255}));
      connect(H.port_b, L6.port_a1) annotation (Line(points={{-360,-46.6},{-352,
              -46.6},{-352,-48},{-346,-48},{-346,-2},{-336,-2}},
                                                               color={0,127,255}));
      connect(G.port_b, L9.port_a1) annotation (Line(points={{-536,-45.8},{-532,
              -45.8},{-532,-46},{-524,-46},{-524,-2},{-518,-2}},
                                                               color={0,127,255}));
      connect(G.port_a, L9.port_b2) annotation (Line(points={{-535.86,-43},{-528,-43},
              {-528,-12},{-518,-12}},      color={0,127,255}));
      connect(P.port_b, L11.port_a1) annotation (Line(points={{-612,-51.8},{
              -606,-51.8},{-606,-52},{-598,-52},{-598,-4},{-594,-4}},
                                                                 color={0,127,255}));
      connect(P.port_a, L11.port_b2) annotation (Line(points={{-611.86,-49},{-604,
              -49},{-604,-14},{-594,-14}}, color={0,127,255}));
      connect(F.port_b, L8.port_a1) annotation (Line(points={{-480,28.2},{-478,
              28.2},{-478,28},{-474,28},{-474,-2},{-456,-2}},
                                                            color={0,127,255}));
      connect(F.port_a, L8.port_b2) annotation (Line(points={{-479.86,31},{-466,31},
              {-466,-12},{-456,-12}},     color={0,127,255}));
      connect(B.port_b, L4.port_a1) annotation (Line(points={{-172,128.2},{-226,
              128.2},{-226,74},{-222,74},{-222,0},{-190,0}},color={0,127,255}));
      connect(B.port_a, L4.port_b2) annotation (Line(points={{-171.86,131},{
              -214,131},{-214,-10},{-190,-10}},
                                          color={0,127,255}));
      connect(E.port_b, L18.port_b1) annotation (Line(points={{-308,108.2},{
              -302,108.2},{-302,106},{-296,106},{-296,130},{-295.8,130}},
                                                                     color={0,127,
              255}));
      connect(E.port_a, L18.port_a2) annotation (Line(points={{-307.86,111},{
              -288,111},{-288,130},{-288.8,130}},
                                             color={0,127,255}));
      connect(D.port_b, L18.port_a1) annotation (Line(points={{-276,174.2},{
              -286,174.2},{-286,174},{-295.8,174},{-295.8,146}},
                                                            color={0,127,255}));
      connect(L18.port_b2, D.port_a) annotation (Line(points={{-288.8,146},{-288,146},
              {-288,180},{-288,179},{-276.14,177}},      color={0,127,255}));
      connect(L13.port_b2, K.port_a) annotation (Line(points={{-26.8,48},{-28,
              48},{-28,83},{-18,83},{1.88,83}},
                                            color={0,127,255}));
      connect(J.port_a, K.port_a) annotation (Line(points={{-45.86,83},{1.88,83}},
                         color={0,127,255}));
      connect(K.port_b, J.port_b) annotation (Line(points={{2,80.2},{-8,80.2},{
              -16,80.2},{-26.8,80.2},{-46,80.2}},                        color={0,
              127,255}));
      connect(O.port_b, L12.port_a1) annotation (Line(points={{-668,14.2},{-672,
              14.2},{-672,14},{-676,14},{-676,-4},{-644,-4}}, color={0,127,255}));
      connect(O.port_a, L12.port_b2) annotation (Line(points={{-668.14,17},{-680,
              17},{-680,-14},{-644,-14}}, color={0,127,255}));
      connect(A.port_a, L15.port_b2) annotation (Line(points={{-237.86,59},{
              -284,59},{-284,56},{-283.2,56},{-283.2,40}}, color={0,127,255}));
      connect(A.port_b, L15.port_a1) annotation (Line(points={{-238,56.2},{-270,
              56.2},{-270,50},{-291.2,50},{-291.2,40}}, color={0,127,255}));
      connect(MassFlowSource.ports[1], L1_2.port_a2) annotation (Line(points={{
              82,-12},{58,-12},{58,-10},{32,-10}}, color={0,127,255}));
      connect(BCTVB.yR,read. u) annotation (Line(
          points={{-853.2,-252},{-788,-252},{-788,230},{-748.4,230}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(write.y,BCTVB. uR) annotation (Line(
          points={{74.4,228},{556,228},{556,-252},{-1101.6,-252}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(A.T_Supply, write.u2[1]) annotation (Line(points={{-254.59,62.71},
              {-170.295,62.71},{-170.295,177.6},{-72.8,177.6}}, color={0,0,127}));
      connect(B.T_Supply, write.u2[2]) annotation (Line(points={{-188.59,134.71},
              {-131.295,134.71},{-131.295,179.2},{-72.8,179.2}}, color={0,0,127}));
      connect(read.y2[5], C.T_re) annotation (Line(points={{-628.8,195.063},{
              -240,195.063},{-240,-53.93},{-231.99,-53.93}}, color={0,0,127}));
      connect(read.y2[6], C.m_flow) annotation (Line(points={{-628.8,195.387},{
              -244,195.387},{-244,-33.14},{-223.8,-33.14}}, color={0,0,127}));
      connect(C.T_Supply, write.u2[3]) annotation (Line(points={{-230.59,-43.29},
              {-250,-43.29},{-250,-66},{-200,-66},{-200,180.8},{-72.8,180.8}},
            color={0,0,127}));
      connect(read.y2[7], D.T_re) annotation (Line(points={{-628.8,195.713},{
              -232,195.713},{-232,168},{-258.01,168},{-258.01,170.07}}, color={
              0,0,127}));
      connect(D.m_flow, read.y2[8]) annotation (Line(points={{-266.2,190.86},{
              -244,190.86},{-244,196.037},{-628.8,196.037}}, color={0,0,127}));
      connect(D.T_Supply, write.u2[4]) annotation (Line(points={{-259.41,180.71},
              {-166.705,180.71},{-166.705,182.4},{-72.8,182.4}}, color={0,0,127}));
      connect(E.T_re, read.y2[9]) annotation (Line(points={{-325.99,104.07},{
              -492,104.07},{-492,196.363},{-628.8,196.363}}, color={0,0,127}));
      connect(E.m_flow, read.y2[10]) annotation (Line(points={{-317.8,124.86},{
              -476,124.86},{-476,196.688},{-628.8,196.688}}, color={0,0,127}));
      connect(E.T_Supply, write.u2[5]) annotation (Line(points={{-324.59,114.71},
              {-328,114.71},{-328,184},{-72.8,184}}, color={0,0,127}));
      connect(F.T_re, read.y2[11]) annotation (Line(points={{-497.99,24.07},{
              -518,24.07},{-518,197.012},{-628.8,197.012}}, color={0,0,127}));
      connect(F.m_flow, read.y2[12]) annotation (Line(points={{-489.8,44.86},{
              -526,44.86},{-526,197.338},{-628.8,197.338}}, color={0,0,127}));
      connect(F.T_Supply, write.u2[6]) annotation (Line(points={{-496.59,34.71},
              {-498,34.71},{-498,185.6},{-72.8,185.6}}, color={0,0,127}));
      connect(G.T_re, read.y2[13]) annotation (Line(points={{-553.99,-49.93},{
              -574,-49.93},{-574,-62},{-514,-62},{-514,197.662},{-628.8,197.662}},
            color={0,0,127}));
      connect(G.m_flow, read.y2[14]) annotation (Line(points={{-545.8,-29.14},{
              -578,-29.14},{-578,-66},{-508,-66},{-508,197.988},{-628.8,197.988}},
            color={0,0,127}));
      connect(G.T_Supply, write.u2[7]) annotation (Line(points={{-552.59,-39.29},
              {-572,-39.29},{-572,-70},{-512,-70},{-512,187.2},{-72.8,187.2}},
            color={0,0,127}));
      connect(H.T_re, read.y2[15]) annotation (Line(points={{-383.13,-51.91},{
              -396,-51.91},{-396,198.313},{-628.8,198.313}}, color={0,0,127}));
      connect(H.m_flow, read.y2[16]) annotation (Line(points={{-372.6,-25.18},{
              -390,-25.18},{-390,198.637},{-628.8,198.637}}, color={0,0,127}));
      connect(H.T_Supply, write.u2[8]) annotation (Line(points={{-381.33,-38.23},
              {-384,-38.23},{-384,188.8},{-72.8,188.8}}, color={0,0,127}));
      connect(I.T_re, read.y2[17]) annotation (Line(points={{-451.99,-79.93},{
              -464,-79.93},{-464,198.963},{-628.8,198.963}}, color={0,0,127}));
      connect(I.m_flow, read.y2[18]) annotation (Line(points={{-443.8,-59.14},{
              -470,-59.14},{-470,199.287},{-628.8,199.287}}, color={0,0,127}));
      connect(I.T_Supply, write.u2[9]) annotation (Line(points={{-450.59,-69.29},
              {-482,-69.29},{-482,190.4},{-72.8,190.4}}, color={0,0,127}));
      connect(J.T_re, read.y2[19]) annotation (Line(points={{-63.99,76.07},{
              -108,76.07},{-108,88},{-554,88},{-554,199.613},{-628.8,199.613}},
            color={0,0,127}));
      connect(J.m_flow, read.y2[20]) annotation (Line(points={{-55.8,96.86},{
              -546,96.86},{-546,199.938},{-628.8,199.938}}, color={0,0,127}));
      connect(K.T_re, read.y2[21]) annotation (Line(points={{17.42,76.07},{40,
              76.07},{40,114},{-564,114},{-564,200.262},{-628.8,200.262}},
            color={0,0,127}));
      connect(K.m_flow, read.y2[22]) annotation (Line(points={{10.4,96.86},{34,
              96.86},{34,118},{-568,118},{-568,200.588},{-628.8,200.588}},
            color={0,0,127}));
      connect(L.T_re, read.y2[23]) annotation (Line(points={{-133.99,-45.93},{
              -184,-45.93},{-184,198},{-628.8,198},{-628.8,200.912}}, color={0,
              0,127}));
      connect(L.m_flow, read.y2[24]) annotation (Line(points={{-125.8,-25.14},{
              -186,-25.14},{-186,201.238},{-628.8,201.238}}, color={0,0,127}));
      connect(L.T_Supply, write.u2[12]) annotation (Line(points={{-132.59,
              -35.29},{-144,-35.29},{-144,195.2},{-72.8,195.2}}, color={0,0,127}));
      connect(J.T_Supply, write.u2[10]) annotation (Line(points={{-62.59,86.71},
              {-108,86.71},{-108,192},{-72.8,192}}, color={0,0,127}));
      connect(K.T_Supply, write.u2[11]) annotation (Line(points={{16.22,86.71},
              {20,86.71},{20,122},{-102,122},{-102,193.6},{-72.8,193.6}}, color=
             {0,0,127}));
      connect(M.T_re, read.y2[25]) annotation (Line(points={{-599.99,94.07},{
              -610,94.07},{-610,201.563},{-628.8,201.563}}, color={0,0,127}));
      connect(M.m_flow, read.y2[26]) annotation (Line(points={{-591.8,114.86},{
              -600,114.86},{-600,201.887},{-628.8,201.887}}, color={0,0,127}));
      connect(M.T_Supply, write.u2[13]) annotation (Line(points={{-598.59,
              104.71},{-602,104.71},{-602,196.8},{-72.8,196.8}}, color={0,0,127}));
      connect(N.T_re, read.y2[27]) annotation (Line(points={{-607.99,42.07},{
              -612,42.07},{-612,202.213},{-628.8,202.213}}, color={0,0,127}));
      connect(N.m_flow, read.y2[28]) annotation (Line(points={{-599.8,62.86},{
              -616,62.86},{-616,202.537},{-628.8,202.537}}, color={0,0,127}));
      connect(N.T_Supply, write.u2[14]) annotation (Line(points={{-606.59,52.71},
              {-608,52.71},{-608,198.4},{-72.8,198.4}}, color={0,0,127}));
      connect(O.T_re, read.y2[29]) annotation (Line(points={{-650.01,10.07},{
              -600,10.07},{-600,202.863},{-628.8,202.863}}, color={0,0,127}));
      connect(O.m_flow, read.y2[30]) annotation (Line(points={{-658.2,30.86},{
              -608,30.86},{-608,203.188},{-628.8,203.188}}, color={0,0,127}));
      connect(O.T_Supply, write.u2[15]) annotation (Line(points={{-651.41,20.71},
              {-528,20.71},{-528,200},{-72.8,200}}, color={0,0,127}));
      connect(P.T_re, read.y2[31]) annotation (Line(points={{-629.99,-55.93},{
              -632,-55.93},{-632,170},{-584,170},{-584,203.512},{-628.8,203.512}},
            color={0,0,127}));
      connect(P.m_flow, read.y2[32]) annotation (Line(points={{-621.8,-35.14},{
              -632,-35.14},{-632,168},{-572,168},{-572,203.838},{-628.8,203.838}},
            color={0,0,127}));
      connect(P.T_Supply, write.u2[16]) annotation (Line(points={{-628.59,
              -45.29},{-594,-45.29},{-594,201.6},{-72.8,201.6}}, color={0,0,127}));
      connect(read.y1[1], MassFlowSource.T_in) annotation (Line(points={{-628.8,
              258.6},{-156,258.6},{-156,50},{116,50},{116,-8},{104,-8}}, color=
              {0,0,127}));
      connect(MassFlowSource.m_flow_in, read.y1[2]) annotation (Line(points={{
              102,-4},{120,-4},{120,-2},{130,-2},{130,60},{-134,60},{-134,263.8},
              {-628.8,263.8}}, color={0,0,127}));
      connect(MassFlowSens.port_a, L1_2.port_b1) annotation (Line(points={{42,16},
              {38,16},{38,0},{32,0}},     color={0,127,255}));
      connect(MassFlowSens.port_b, TempSens.port_a) annotation (Line(points={{62,16},
              {108,16},{108,18},{142,18}},        color={0,127,255}));
      connect(TempSens.port_b, Sink.ports[1]) annotation (Line(points={{162,18},
              {200,18},{200,12}}, color={0,127,255}));
      connect(MassFlowSens.m_flow, write.u1[1]) annotation (Line(points={{52,27},
              {60,27},{60,142},{-120,142},{-120,268},{-72.8,268},{-72.8,260}},
            color={0,0,127}));
      connect(TempSens.T, write.u1[2]) annotation (Line(points={{152,29},{150,
              29},{150,158},{-112,158},{-112,272.8},{-72.8,272.8}}, color={0,0,
              127}));
      connect(N.T_Supply, write.u2[14]) annotation (Line(points={{-606.59,52.71},
              {-610,52.71},{-610,198.4},{-72.8,198.4}}, color={0,0,127}));
      connect(L13.port_a1, J.port_b) annotation (Line(points={{-18.8,48},{-20,
              48},{-20,80},{-26.8,80.2},{-46,80.2}}, color={0,127,255}));
      connect(K.port_b, L13.port_a1) annotation (Line(points={{2,80.2},{-10,
              80.2},{-10,80},{-18.8,80},{-18.8,48}}, color={0,127,255}));
      connect(L13.port_b2, J.port_a) annotation (Line(points={{-26.8,48},{-26.8,
              83},{-45.86,83}}, color={0,127,255}));
      connect(A.T_re, read.y2[1]) annotation (Line(points={{-255.99,52.07},{
              -354,52.07},{-354,193.762},{-628.8,193.762}}, color={0,0,127}));
      connect(A.m_flow, read.y2[2]) annotation (Line(points={{-247.8,72.86},{
              -404,72.86},{-404,194.088},{-628.8,194.088}}, color={0,0,127}));
      connect(B.T_re, read.y2[3]) annotation (Line(points={{-189.99,124.07},{
              -230,124.07},{-230,194.412},{-628.8,194.412}}, color={0,0,127}));
      connect(B.m_flow, read.y2[4]) annotation (Line(points={{-181.8,144.86},{
              -230,144.86},{-230,194.738},{-628.8,194.738}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-1180,
                -520},{640,340}})), Diagram(coordinateSystem(preserveAspectRatio=
                false, extent={{-1180,-520},{640,340}})));
    end GrazNetworkModel2;
  end Experiments;
  annotation (uses(
      Modelica(version="3.2.2"),
      Annex60(version="0.1"),
      Buildings(version="3.0.0"),
      DHNLibrary(version="3")));
end BCVTB_Graz_Test;
