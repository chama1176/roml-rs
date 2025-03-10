#[allow(unused_imports)]

/// Robot Kinematics and Dynamics
pub struct IK4dTriangle<T> {
    pub a: T,
    pub b: T,
}

impl<T: na::RealField> IK4dTriangle<T> {
    pub fn new() -> Self {
        Self {
            a: T::zero(),
            b: T::zero(),
        }
    }

    pub fn solve_by_ref_point(&self, vc: &na::Vector3<T>, ref_point: &na::Vector3<T>) -> [T; 4] {
        let mut ans = [T::zero(), T::zero(), T::zero(), T::zero()];

        let c_squared =
            vc.x.clone() * vc.x.clone() + vc.y.clone() * vc.y.clone() + vc.z.clone() * vc.z.clone();
        let c = c_squared.sqrt();
        // cosine formula
        ans[3] = T::pi()
            - ((self.a.clone() * self.a.clone() + self.b.clone() * self.b.clone()
                - c.clone() * c.clone())
                / (T::from_f32(2.0).unwrap() * self.a.clone() * self.b.clone()))
            .acos();

        // angle between vc and va
        let b_theta = ((self.a.clone() * self.a.clone() + c.clone() * c.clone()
            - self.b.clone() * self.b.clone())
            / (T::from_f32(2.0).unwrap() * self.a.clone() * c.clone()))
        .acos();

        // rotation from axis vertial to ref point vector and vc
        let rot_c_ref = na::UnitQuaternion::from_axis_angle(
            &na::Unit::new_normalize(vc.cross(&ref_point)),
            b_theta.clone(),
        );
        let va = rot_c_ref.transform_vector(&vc) * self.a.clone() / c.clone();
        let vb = vc.clone() - va.clone();

        // // rotation along with vector c
        // let rot_ref_theta = na::UnitQuaternion::from_axis_angle(
        //     &na::Unit::new_normalize(vc.clone()),
        //     ref_theta.clone(),
        // );
        // let va = rot_ref_theta.transform_vector(&va_dot);
        // let vb = rot_ref_theta.transform_vector(&vb_dot);
        ans[0] = va.y.clone().atan2(va.x.clone());
        ans[1] = (va.z.clone() / self.a.clone()).acos(); // 精度悪いめ

        let rot_ans0 = na::UnitQuaternion::from_axis_angle(&na::Vector3::z_axis(), ans[0].clone());
        let rot_ans3 = na::UnitQuaternion::from_axis_angle(
            &na::Unit::new_normalize(rot_ans0.transform_vector(&na::Vector3::y_axis())),
            ans[3].clone(),
        );

        let vb_dot_dot = rot_ans3.transform_vector(&va) * self.b.clone() / self.a.clone();

        // vector a component of vector b
        let va_tmp = va.clone() * va.dot(&vb) / self.a.clone() / self.a.clone();
        // let va_tmp2 = va.clone() * va.dot(&vb_dot_dot) / self.a.clone() / self.a.clone();
        // assert_eq!(va_tmp, va_tmp2);

        let vb_tmp = vb_dot_dot.clone() - va_tmp.clone();
        let vb_tmp2 = vb.clone() - va_tmp.clone();

        let vb_tmp_cross = vb_tmp.cross(&vb_tmp2);
        ans[2] = va.clone().dot(&vb_tmp_cross).signum()
            * vb_tmp_cross.dot(&vb_tmp_cross).sqrt().atan2(vb_tmp.dot(&vb_tmp2));

        ans
        // 👺原点と方向の定義も必要
    }
    pub fn solve(&self, vc: &na::Vector3<T>, ref_theta: &T) -> [T; 4] {
        
        let rot_c_ref = na::UnitQuaternion::from_axis_angle(
            &na::Unit::new_normalize(vc.clone()),
            ref_theta.clone(),
        );
        let ref_point = rot_c_ref.transform_vector(&na::Vector3::z_axis());

        self.solve_by_ref_point(vc, &ref_point)

    }
}


pub struct IK4dExtendedTriangle<T> {
    pub a: T,
    pub b: T,
    ik4dtriangle :IK4dTriangle<T>,
}

impl<T: na::RealField> IK4dExtendedTriangle<T> {
    pub fn new() -> Self {
        Self {
            a: T::zero(),
            b: T::zero(),
            ik4dtriangle: IK4dTriangle::new(),
        }
    }

}

#[cfg(test)]
mod test_ik {
    use crate::ik::IK4dTriangle;
    use approx::assert_relative_eq;

    trait Deg {
        fn deg2rad(&self) -> Self;
    }

    impl Deg for f32 {
        fn deg2rad(&self) -> Self {
            self * core::f32::consts::PI / 180.
        }
    }

    #[test]
    fn ik_4dof_equilateral_triangle() {
        let mut ik = IK4dTriangle::<f32>::new();

        ik.a = 3.0;
        ik.b = 3.0;
        let ans = ik.solve(&na::Vector3::new(3.0, 0.0, 0.0), &0.0);
        assert_relative_eq!(ans[0], 0.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[1], core::f32::consts::PI / 6.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[2], 0.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[3], core::f32::consts::PI * 2.0 / 3.0, epsilon = 1.0e-6);
        ik.a = 3.0;
        ik.b = 3.0;
        let ans = ik.solve(&na::Vector3::new(
            3.0 / (2.0 as f32).sqrt(),
            3.0 / (2.0 as f32).sqrt(),
            0.0,
        ), &0.0);
        assert_relative_eq!(ans[0], core::f32::consts::PI / 4.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[1], core::f32::consts::PI / 6.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[2], 0.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[3], core::f32::consts::PI * 2.0 / 3.0, epsilon = 1.0e-6);
        ik.a = 3.0;
        ik.b = 3.0;
        let ans = ik.solve(&na::Vector3::new(3.0, 0.0, 0.0), &(core::f32::consts::PI / 2.0));
        assert_relative_eq!(ans[0], -core::f32::consts::PI / 3.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[1], core::f32::consts::PI / 2.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[2], core::f32::consts::PI / 2.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[3], core::f32::consts::PI * 2.0 / 3.0, epsilon = 1.0e-6);
        ik.a = 3.0;
        ik.b = 3.0;
        let ans = ik.solve(&na::Vector3::new(3.0, 0.0, 0.0), &(-core::f32::consts::PI / 2.0));
        assert_relative_eq!(ans[0], core::f32::consts::PI / 3.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[1], core::f32::consts::PI / 2.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[2], -core::f32::consts::PI / 2.0, epsilon = 1.0e-6);
        assert_relative_eq!(ans[3], core::f32::consts::PI * 2.0 / 3.0, epsilon = 1.0e-6);
    }

    #[test]
    fn ik_4dof_triangle() {
        let mut ik = IK4dTriangle::<f32>::new();
        ik.a = 1.0;
        ik.b = (3.0 as f32).sqrt();
        let ans = ik.solve(&na::Vector3::new(2.0, 0.0, 0.0), &(core::f32::consts::PI / 4.));
        assert_relative_eq!(ans[0], -50.76847952.deg2rad(), epsilon = 1.0e-6);
        assert_relative_eq!(ans[1], 52.23875609.deg2rad(), epsilon = 1.0e-6);
        assert_relative_eq!(ans[2], 63.43494882.deg2rad(), epsilon = 1.0e-6);
        assert_relative_eq!(ans[3], core::f32::consts::PI / 2.0, epsilon = 1.0e-6);
    }
}
