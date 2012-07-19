/** demos.h
	Copyright (C) 2012 Jackson Lee

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/
#ifndef FLUID_SYSTEM_DEMOS_H_INCLUDED
#define FLUID_SYSTEM_DEMOS_H_INCLUDED

#include "btBulletDynamicsCommon.h"
#include "Fluids/FluidSph.h"
#include "Fluids/FluidWorld.h"
#include "Fluids/SPHInterface.h"

#include "BulletCollision/CollisionShapes/btTriangleMesh.h"

class FluidSystemDemo
{
protected:
	static btRigidBody* createRigidBody(const btTransform &transform, btScalar mass, btCollisionShape *shape)
	{
		//Rigid bodies are dynamic if and only if mass is non zero, otherwise static
		btVector3 localInertia(0,0,0);
		if(mass != 0.f) shape->calculateLocalInertia(mass, localInertia);

		btDefaultMotionState *motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
		btRigidBody *result = new btRigidBody(rbInfo);
		
		return result;
	}

	btAlignedObjectArray<btRigidBody*> m_rigidBodies;

public:
	virtual ~FluidSystemDemo() {}
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes) {}
	virtual void deactivate() 
	{
		for(int i = 0; i < m_rigidBodies.size(); ++i)
		{
			if(m_rigidBodies[i])
			{
				if( m_rigidBodies[i]->getMotionState() ) delete m_rigidBodies[i]->getMotionState();
				delete m_rigidBodies[i];
			}
		}
		
		m_rigidBodies.clear();
	}
	
	virtual void addToWorld(btDynamicsWorld *world) { for(int i = 0; i < m_rigidBodies.size(); ++i)world->addRigidBody(m_rigidBodies[i]); }
	virtual void removeFromWorld(btDynamicsWorld *world) { for(int i = 0; i < m_rigidBodies.size(); ++i)world->removeRigidBody(m_rigidBodies[i]); }
	
	//	rename as 'stepSimulation'?
	virtual void update(const FluidWorld &FW, FluidSph *fluid) {}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles) = 0;
};
class Demo_DamBreak : public FluidSystemDemo
{
public:
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 20.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
		
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-INIT_BOUND, 10.f, 0.f);
		btVector3 initMax(INIT_BOUND, 55.f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	}
};
class Demo_Drop : public FluidSystemDemo
{
public:
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 30.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
	
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		btVector3 initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	}
};
class Demo_EmitterAndAbsorber : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 0.0;

		btCollisionShape *wallShape = new btBoxShape( btVector3(5.0, 15.0, 50.0) );
		collisionShapes->push_back(wallShape);
		
		btTransform transform;
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(25.0, 10.0, 0.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(-25.0, 10.0, 0.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(0.0, 10.0, 25.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(55.0, 10.0, -25.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(-55.0, 10.0, -25.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
	}
	
	virtual void update(const FluidWorld &FW, FluidSph *fluid)
	{
		FluidEmitter emitter;
		emitter.m_pitch = -90.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(10.f, 10.f, 10.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 20;
			emitter.emit( fluid, NUM_PARTICLES_EMITTED, fluid->getEmitterSpacing(FW.getGlobalParameters()) );
		}
		
		//
		FluidAbsorber absorber;
		absorber.m_min.setValue(-100.0, -20.0, -75.0);
		absorber.m_max.setValue(100.0, 20.0, -50.0);
		absorber.absorb(fluid);
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 120.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
	
		//const btScalar INIT_BOUND = 20.0f;
		//btVector3 initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		//btVector3 initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		//FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	}
};

class Demo_Levee : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 0.0;
		
		btTransform transform;
		
		//Sloped ground
		btCollisionShape* squareBoxShape = new btBoxShape( btVector3(50.0, 5.0, 50.0) );
		collisionShapes->push_back(squareBoxShape);
		
		transform = btTransform( btQuaternion(0.0, 0.0, 0.1), btVector3(0.0, 0.0, 0.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, squareBoxShape) );
		
		//Walls
		btCollisionShape *wallShape = new btBoxShape( btVector3(1.0, 5.0, 30.0) );
		//btCollisionShape *wallShape = new btBoxShape( btVector3(0.1, 10.0, 30.0) );	//for CCD_TEST
		collisionShapes->push_back(wallShape);
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(0.0, 5.0, -35.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(0.0, 5.0, 35.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
	}
	
	virtual void update(const FluidWorld &FW, FluidSph *fluid)
	{
		FluidEmitter emitter;
		emitter.m_pitch = 225.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(35.f, 20.f, -20.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 3;
			emitter.emit( fluid, NUM_PARTICLES_EMITTED, fluid->getEmitterSpacing(FW.getGlobalParameters()) );
		}
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 40.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
		
		const btScalar INIT_BOUND = 30.0f;
		btVector3 initMin(20.0, 20.0f, -INIT_BOUND);
		btVector3 initMax(INIT_BOUND, 70.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	}
};


class Demo_Drain : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 0.0;
		
		const btScalar TILE_EXTENT = 10.0;
		const btScalar TILE_POSITION = TILE_EXTENT * 2.0;
		const btScalar HEIGHT = 20.0;
		
		btTransform transform;
		
		btCollisionShape* tileShape = new btBoxShape( btVector3(TILE_EXTENT, 3.0, TILE_EXTENT) );
		collisionShapes->push_back(tileShape);

		//Create 8 tiles
		transform = btTransform( btQuaternion::getIdentity(), btVector3(TILE_POSITION, HEIGHT, 0.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(-TILE_POSITION, HEIGHT, 0.0) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(0.0, HEIGHT, TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(0.0, HEIGHT, -TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(TILE_POSITION, HEIGHT, TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(TILE_POSITION, HEIGHT, -TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(-TILE_POSITION, HEIGHT, TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(-TILE_POSITION, HEIGHT, -TILE_POSITION) );
		m_rigidBodies.push_back( createRigidBody(transform, MASS, tileShape) );
	}
	
	virtual void update(const FluidWorld &FW, FluidSph *fluid)
	{
		FluidEmitter emitter;
		emitter.m_pitch = 225.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(25.f, 30.f, -15.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 3;
			emitter.emit( fluid, NUM_PARTICLES_EMITTED, fluid->getEmitterSpacing(FW.getGlobalParameters()) );
		}
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 30.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*3.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
	
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-INIT_BOUND, 30.0f, -INIT_BOUND);
		btVector3 initMax(INIT_BOUND, 70.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	}
};

class Demo_DynamicBox : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 1.f;
	
		btCollisionShape* largeBoxShape = new btBoxShape( btVector3(10.0, 10.0, 10.0) );
		collisionShapes->push_back(largeBoxShape);
		
		btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 10.0, 0.0) );
		m_rigidBodies.push_back( createRigidBody(startTransform, MASS, largeBoxShape) );
	}
	
	virtual void update(const FluidWorld &FW, FluidSph *fluid)
	{
		FluidEmitter emitter;
		emitter.m_pitch = -90.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(10.f, 10.f, 10.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 5;
			emitter.emit( fluid, NUM_PARTICLES_EMITTED, fluid->getEmitterSpacing(FW.getGlobalParameters()) );
		}
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 30.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
		
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		btVector3 initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
	
		if( m_rigidBodies.size() && m_rigidBodies[0] )
		{
			btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 10.0, 0.0) );
			m_rigidBodies[0]->setCenterOfMassTransform(startTransform);
		}
	}
};

	//	next: add a bucket composed of constrained btBoxShapes
class Demo_HollowBox : public FluidSystemDemo
{
	btTriangleMesh *m_hollowBoxMesh;

public:
	Demo_HollowBox() : m_hollowBoxMesh(0) {}

	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		if(!m_hollowBoxMesh)
		{
			const btScalar BUCKET_SIZE = 8.0;
			m_hollowBoxMesh = new btTriangleMesh(true, true);	
			
			//RL: Right/Left == +X/-X
			//UD: Up/Down == +Y/-Y
			//BF: Backwards/Forwards == +Z/-Z
			btVector3 RUB(BUCKET_SIZE, BUCKET_SIZE, BUCKET_SIZE);
			btVector3 RUF(BUCKET_SIZE, BUCKET_SIZE, -BUCKET_SIZE);
			btVector3 RDB(BUCKET_SIZE, -BUCKET_SIZE, BUCKET_SIZE);
			btVector3 RDF(BUCKET_SIZE, -BUCKET_SIZE, -BUCKET_SIZE);
			btVector3 LUB(-BUCKET_SIZE, BUCKET_SIZE, BUCKET_SIZE);
			btVector3 LUF(-BUCKET_SIZE, BUCKET_SIZE, -BUCKET_SIZE);
			btVector3 LDB(-BUCKET_SIZE, -BUCKET_SIZE, BUCKET_SIZE);
			btVector3 LDF(-BUCKET_SIZE, -BUCKET_SIZE, -BUCKET_SIZE);
			
			//Face +X
			m_hollowBoxMesh->addTriangle(RDB, RUF, RUB, true);
			m_hollowBoxMesh->addTriangle(RDB, RUF, RDF, true);
			
			//Face -X
			m_hollowBoxMesh->addTriangle(LDB, LUF, LUB, true);
			m_hollowBoxMesh->addTriangle(LDB, LUF, LDF, true);
			
			//Face +Y
			m_hollowBoxMesh->addTriangle(RUB, LUF, RUF, true);
			m_hollowBoxMesh->addTriangle(RUB, LUF, LUB, true);
			
			//Face -Y
			m_hollowBoxMesh->addTriangle(RDB, LDF, RDF, true);
			m_hollowBoxMesh->addTriangle(RDB, LDF, LDB, true);
			
			//Face +Z
			m_hollowBoxMesh->addTriangle(RUB, LDB, LUB, true);
			m_hollowBoxMesh->addTriangle(RUB, LDB, RDB, true);
			
			//Face -Z
			m_hollowBoxMesh->addTriangle(RUF, LDF, LUF, true);
			m_hollowBoxMesh->addTriangle(RUF, LDF, RDF, true);
		}
		
		btBvhTriangleMeshShape *bucketShape = new btBvhTriangleMeshShape(m_hollowBoxMesh, true);
		collisionShapes->push_back(bucketShape);
		
		const btScalar MASS = 1.0;
		btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 8.0, 0.0) );
		m_rigidBodies.push_back( createRigidBody(startTransform, MASS, bucketShape) );
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 30.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
		
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		btVector3 initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
		
		if( m_rigidBodies.size() && m_rigidBodies[0] )
		{
			btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 8.0, 0.0) );
			m_rigidBodies[0]->setCenterOfMassTransform(startTransform);
		}
	}
	
	virtual void deactivate()
	{
		FluidSystemDemo::deactivate();
		
		if(m_hollowBoxMesh)
		{
			delete m_hollowBoxMesh;
			m_hollowBoxMesh = 0;
		}
	}
};
	
class Demo_FluidToRigidBody : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		btCollisionShape *wallShape = new btBoxShape( btVector3(5.0, 20.0, 50.0) );
		collisionShapes->push_back(wallShape);
		{
			const btScalar MASS = 0.0;

			btTransform transform;
			
			transform = btTransform( btQuaternion::getIdentity(), btVector3(27.0, 10.0, 0.0) );
			m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
			
			transform = btTransform( btQuaternion::getIdentity(), btVector3(-27.0, 10.0, 0.0) );
			m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
			
			transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(0.0, 10.0, 50.0) );
			m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
			
			transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(0.0, 10.0, -50.0) );
			m_rigidBodies.push_back( createRigidBody(transform, MASS, wallShape) );
		}
		
		btCollisionShape* boxShape = new btBoxShape( btVector3(3.0, 3.0, 3.0) );
		collisionShapes->push_back(boxShape);
		{
			const btScalar MASS = 0.01f;
		
			btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 3.0, 0.0) );
			m_rigidBodies.push_back( createRigidBody(startTransform, MASS, boxShape) );
			btTransform startTransform2( btQuaternion::getIdentity(), btVector3(0.0, 9.0, 0.0) );
			m_rigidBodies.push_back( createRigidBody(startTransform2, MASS, boxShape) );
			btTransform startTransform3( btQuaternion::getIdentity(), btVector3(0.0, 15.0, 0.0) );
			m_rigidBodies.push_back( createRigidBody(startTransform3, MASS, boxShape) );
		}
	}
	
	virtual void reset(const FluidWorld &FW, FluidSph *fluid, int maxFluidParticles)
	{
		const btScalar VOL_BOUND = 50.0f;
		btVector3 volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		btVector3 volMax(VOL_BOUND, 80.0f, VOL_BOUND);
		
		fluid->initialize(FW.getGlobalParameters(), maxFluidParticles, volMin, volMax);
	
		const btScalar INIT_BOUND = 20.0f;
		btVector3 initMin(-20.0f, 20.0f, 20.0f);
		btVector3 initMax(20.0f, 60.0f, 40.0f);
		FluidEmitter::addVolume( fluid, initMin, initMax, fluid->getEmitterSpacing(FW.getGlobalParameters()) * 0.87 );
		
		if( m_rigidBodies.size() )
		{
			if(m_rigidBodies[4])
			{
				btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 3.0, 0.0) );
				m_rigidBodies[4]->setCenterOfMassTransform(startTransform);
				m_rigidBodies[4]->setLinearVelocity( btVector3(0,0,0) );
				m_rigidBodies[4]->setAngularVelocity( btVector3(0,0,0) );
			}
			if(m_rigidBodies[5])
			{
				btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 9.0, 0.0) );
				m_rigidBodies[5]->setCenterOfMassTransform(startTransform);
				m_rigidBodies[5]->setLinearVelocity( btVector3(0,0,0) );
				m_rigidBodies[5]->setAngularVelocity( btVector3(0,0,0) );
			}
			if(m_rigidBodies[6])
			{
				btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 15.0, 0.0) );
				m_rigidBodies[6]->setCenterOfMassTransform(startTransform);
				m_rigidBodies[6]->setLinearVelocity( btVector3(0,0,0) );
				m_rigidBodies[6]->setAngularVelocity( btVector3(0,0,0) );
			}
		}
	}
};

#endif